import pandas as pd
import os
import glob
import numpy as np
import hashlib
import subprocess
from typing import List, Dict

# ────────────────────────────────────────────────────────────────
# Configuration
# ────────────────────────────────────────────────────────────────

CZID_DIR = 'czid'
SEQTOID_DIR = 'seqtoid'

EXPECTED_SAMPLES = [
    'SRR18291896',
    'SRR15049352',
    'SRR12048509'
]

DIFF_SYMBOLS = {
    'equivalent':   '✅ <0.005',
    'warning':      '⚠️ 0.005–0.05',
    'significant':  '❌ >0.05',
    'no diffs':     'identical',
    'identical':    'T',
    'differ':       'F',
    'empty':        'EMPTY (0 contigs)'
}

# ────────────────────────────────────────────────────────────────
# Helpers
# ────────────────────────────────────────────────────────────────

def file_sha256(filepath):
    """Compute SHA-256 hash of a file in chunks."""
    sha256 = hashlib.sha256()
    with open(filepath, 'rb') as f:
        while chunk := f.read(8192 * 1024):
            sha256.update(chunk)
    return sha256.hexdigest()

def get_fasta_total_length(fasta_path):
    """
    Calculate total base pairs in a FASTA file by summing lengths of sequence lines.
    Ignores headers and whitespace/newlines.
    """
    total = 0
    try:
        with open(fasta_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('>'):
                    total += len(line)
        return total
    except Exception as e:
        print(f"  Warning: Could not read {fasta_path} for length: {e}")
        return 0


def numeric_diff(a: np.ndarray, b: np.ndarray, atol: float = 0.005) -> str:
    """Categorize worst difference."""
    if len(a) == 0 or len(b) == 0:
        return "empty"
    diff = np.abs(a - b)
    worst = diff.max()
    if np.isnan(worst):
        return "NaN"
    if worst <= 0.005:
        return 'equivalent'
    elif worst <= 0.05:
        return 'warning'
    else:
        return 'significant'


def compare_numeric_dfs(df1: pd.DataFrame, df2: pd.DataFrame, id_cols: List[str]) -> Dict[str, str]:
    num_cols = df1.select_dtypes(include=[np.number]).columns.intersection(df2.columns)
    results = {}
    for col in num_cols:
        vals1 = df1[col].fillna(0).values.astype(float)
        vals2 = df2[col].fillna(0).values.astype(float)
        cat = numeric_diff(vals1, vals2)
        results[col] = cat
    return results


# ────────────────────────────────────────────────────────────────
# Step 1: Sample Metadata
# ────────────────────────────────────────────────────────────────

def compare_metadata():
    file = 'sample_metadata.csv'
    czid_path = os.path.join(CZID_DIR, file)
    seqtoid_path = os.path.join(SEQTOID_DIR, file)

    if not os.path.exists(czid_path) or not os.path.exists(seqtoid_path):
        print(f"Error: {file} missing in one or both dirs.")
        pd.DataFrame({'file': [file], 'identical': ['missing']}).to_csv('step1_metadata.csv', index=False)
        return

    czid_df = pd.read_csv(czid_path)
    seqtoid_df = pd.read_csv(seqtoid_path)

    czid_s = czid_df.sort_values('sample_name').reset_index(drop=True)
    seqtoid_s = seqtoid_df.sort_values('sample_name').reset_index(drop=True)

    identical = czid_s.equals(seqtoid_s)
    print(f"Step 1: sample_metadata.csv → {'identical' if identical else 'differ'}")

    pd.DataFrame({'file': [file], 'identical': ['T' if identical else 'F']}).to_csv('step1_metadata.csv', index=False)

    # Missing / extra samples
    rows = []
    for name, df in [('czid', czid_df), ('seqtoid', seqtoid_df)]:
        actual = set(df['sample_name'].astype(str))
        missing = ', '.join(sorted(set(EXPECTED_SAMPLES) - actual))
        extra = ', '.join(sorted(actual - set(EXPECTED_SAMPLES)))
        if missing or extra:
            rows.append({'directory': name, 'missing': missing, 'extra': extra})
    if rows:
        pd.DataFrame(rows).to_csv('step1_missing_extra.csv', index=False)


# Step 2: Sample Overviews
# (unchanged from your version)

def compare_sample_overviews():
    file_name = 'sample_overviews.csv'
    czid_path = os.path.join(CZID_DIR, file_name)
    seqtoid_path = os.path.join(SEQTOID_DIR, file_name)

    if not os.path.exists(czid_path) or not os.path.exists(seqtoid_path):
        print(f"Error: {file_name} missing in one or both dirs.")
        pd.DataFrame({'file': [file_name], 'strict_identical': ['missing']}).to_csv('step2_overviews.csv', index=False)
        return

    czid_df = pd.read_csv(czid_path)
    seqtoid_df = pd.read_csv(seqtoid_path)

    czid_s = czid_df.sort_values('sample_name').reset_index(drop=True)
    seqtoid_s = seqtoid_df.sort_values('sample_name').reset_index(drop=True)

    strict_identical = czid_s.equals(seqtoid_s)

    results = {'file': file_name, 'strict_identical': 'T' if strict_identical else 'F'}

    if not strict_identical:
        num_results = compare_numeric_dfs(czid_s, seqtoid_s, id_cols=['sample_name'])
        for col, cat in num_results.items():
            results[col] = DIFF_SYMBOLS.get(cat, cat)

    pd.DataFrame([results]).to_csv('step2_overviews.csv', index=False)
    print(f"Step 2: sample_overviews.csv → {'identical' if strict_identical else 'differ'}")


# Step 3: Per-sample Taxon Reports (unchanged)

def compare_taxon_reports():
    print("\n=== Step 3: Sample Taxon Reports Comparison ===")
    rows = []
    missing_czid = []
    missing_seqtoid = []

    for sample in EXPECTED_SAMPLES:
        czid_files = glob.glob(os.path.join(CZID_DIR, f"{sample}_*_taxon_report.csv"))
        seqtoid_files = glob.glob(os.path.join(SEQTOID_DIR, f"{sample}_*_taxon_report.csv"))

        if len(czid_files) != 1 or len(seqtoid_files) != 1:
            if len(czid_files) != 1: missing_czid.append(sample)
            if len(seqtoid_files) != 1: missing_seqtoid.append(sample)
            rows.append({'sample': sample, 'status': 'missing file(s)'})
            continue

        czid_path = czid_files[0]
        seqtoid_path = seqtoid_files[0]

        try:
            czid_df = pd.read_csv(czid_path, dtype=str)
            seqtoid_df = pd.read_csv(seqtoid_path, dtype=str)

            czid_df['tax_id'] = pd.to_numeric(czid_df['tax_id'], errors='coerce')
            seqtoid_df['tax_id'] = pd.to_numeric(seqtoid_df['tax_id'], errors='coerce')

            # Sort by tax_id (stable key)
            czid_s = czid_df.sort_values('tax_id').reset_index(drop=True)
            seqtoid_s = seqtoid_df.sort_values('tax_id').reset_index(drop=True)

            result_row = {'sample': sample}

            # Always run tolerant numeric comparison
            num_results = compare_numeric_dfs(czid_s, seqtoid_s, id_cols=['tax_id'])

            if not num_results:
                result_row['overall'] = 'no numeric columns'
            else:
                # Summarize worst category across all columns
                worst_cat = max(num_results.values(), key=lambda c: ['equivalent', 'warning', 'significant'].index(c) if c in ['equivalent', 'warning', 'significant'] else -1)
                result_row['overall'] = DIFF_SYMBOLS.get(worst_cat, worst_cat)

                # Include per-column results
                for col, cat in num_results.items():
                    result_row[col] = DIFF_SYMBOLS.get(cat, cat)

            # Also report row count difference if present
            if len(czid_s) != len(seqtoid_s):
                result_row['row_count_diff'] = f"czid: {len(czid_s)}, seqtoid: {len(seqtoid_s)}"

            rows.append(result_row)

        except Exception as e:
            rows.append({'sample': sample, 'status': f'error: {str(e)}'})

    # Write main results
    df_results = pd.DataFrame(rows)
    df_results.to_csv('step3_taxon_reports.csv', index=False)

    # Missing files summary
    if missing_czid or missing_seqtoid:
        pd.DataFrame({
            'missing_czid': [', '.join(sorted(missing_czid))],
            'missing_seqtoid': [', '.join(sorted(missing_seqtoid))]
        }).to_csv('step3_missing.csv', index=False)

    print(f"  → Results saved to step3_taxon_reports.csv")


# Step 4: Combined Taxon Results (unchanged)

def compare_combined_taxon_results():
    file_name = 'combined_sample_taxon_results_NT.bpm.csv'
    czid_path = os.path.join(CZID_DIR, file_name)
    seqtoid_path = os.path.join(SEQTOID_DIR, file_name)

    if not os.path.exists(czid_path) or not os.path.exists(seqtoid_path):
        print(f"Error: {file_name} missing in one or both dirs.")
        pd.DataFrame({'file': [file_name], 'strict_identical': ['missing']}).to_csv('step4_combined_taxon_bpm.csv', index=False)
        return

    czid_df = pd.read_csv(czid_path)
    seqtoid_df = pd.read_csv(seqtoid_path)

    sort_col = 'Taxon Name'
    czid_s = czid_df.sort_values(sort_col).reset_index(drop=True)
    seqtoid_s = seqtoid_df.sort_values(sort_col).reset_index(drop=True)

    strict_identical = czid_s.equals(seqtoid_s)

    results = {'file': file_name, 'strict_identical': 'T' if strict_identical else 'F'}

    if not strict_identical:
        num_results = compare_numeric_dfs(czid_s, seqtoid_s, id_cols=[sort_col])
        for col, cat in num_results.items():
            results[col] = DIFF_SYMBOLS.get(cat, cat)

    pd.DataFrame([results]).to_csv('step4_combined_taxon_bpm.csv', index=False)

    print(f"Step 4: combined_sample_taxon_results_NT.bpm.csv → {'identical' if strict_identical else 'differ'}")


# ────────────────────────────────────────────────────────────────
# Step 5: Contig Summary Reports – FIXED for 0-byte / empty files
# ────────────────────────────────────────────────────────────────

def compare_contig_summary_reports():
    print("\n=== Step 5: Contig Summary Reports Comparison ===")
    print("    → Comparing only NT.species_taxid and NR.species_taxid per contig")

    rows = []
    missing_czid = []
    missing_seqtoid = []

    for sample in EXPECTED_SAMPLES:
        czid_files = glob.glob(os.path.join(CZID_DIR, f"{sample}_*_contig_summary_report.csv"))
        seqtoid_files = glob.glob(os.path.join(SEQTOID_DIR, f"{sample}_*_contig_summary_report.csv"))

        if len(czid_files) != 1 or len(seqtoid_files) != 1:
            if len(czid_files) != 1: missing_czid.append(sample)
            if len(seqtoid_files) != 1: missing_seqtoid.append(sample)
            rows.append({'sample': sample, 'status': 'missing file(s)'})
            continue

        czid_path = czid_files[0]
        seqtoid_path = seqtoid_files[0]

        czid_size = os.path.getsize(czid_path)
        seqtoid_size = os.path.getsize(seqtoid_path)

        # Both empty → treated as match (no contigs)
        if czid_size == 0 and seqtoid_size == 0:
            rows.append({'sample': sample, 'status': 'both empty (0 contigs)', 'match': 'T'})
            continue

        # One empty, one not → mismatch
        if czid_size == 0 or seqtoid_size == 0:
            rows.append({
                'sample': sample,
                'status': f'one side empty (czid={czid_size} bytes, seqtoid={seqtoid_size} bytes)',
                'match': 'F'
            })
            continue

        try:
            czid_df = pd.read_csv(czid_path, dtype=str)
            seqtoid_df = pd.read_csv(seqtoid_path, dtype=str)
        except pd.errors.EmptyDataError:
            rows.append({'sample': sample, 'status': 'empty file (read error)', 'match': 'F'})
            continue

        # Only keep the columns we care about
        key_cols = ['contig_name', 'NT.species_taxid', 'NR.species_taxid']
        missing_cols_czid = [c for c in key_cols if c not in czid_df.columns]
        missing_cols_seqtoid = [c for c in key_cols if c not in seqtoid_df.columns]

        if missing_cols_czid or missing_cols_seqtoid:
            rows.append({
                'sample': sample,
                'status': f'missing columns: czid={missing_cols_czid}, seqtoid={missing_cols_seqtoid}',
                'match': 'F'
            })
            continue

        # Select only relevant columns
        czid_df = czid_df[key_cols].copy()
        seqtoid_df = seqtoid_df[key_cols].copy()

        # Clean up taxid columns (handle NaN, convert to str for exact match)
        for df in [czid_df, seqtoid_df]:
            df['NT.species_taxid'] = df['NT.species_taxid'].fillna('NA').astype(str)
            df['NR.species_taxid'] = df['NR.species_taxid'].fillna('NA').astype(str)

        # Try to sort/match by contig_name
        sort_col = 'contig_name'
        if sort_col in czid_df.columns and sort_col in seqtoid_df.columns:
            czid_df[sort_col] = czid_df[sort_col].astype(str).str.strip()
            seqtoid_df[sort_col] = seqtoid_df[sort_col].astype(str).str.strip()

            czid_s = czid_df.sort_values(sort_col).reset_index(drop=True)
            seqtoid_s = seqtoid_df.sort_values(sort_col).reset_index(drop=True)
        else:
            # fallback: just reset index (order as in file)
            czid_s = czid_df.reset_index(drop=True)
            seqtoid_s = seqtoid_df.reset_index(drop=True)

        # Compare row counts first
        if len(czid_s) != len(seqtoid_s):
            rows.append({
                'sample': sample,
                'status': f'different contig count (czid: {len(czid_s)}, seqtoid: {len(seqtoid_s)})',
                'match': 'F'
            })
            continue

        # Exact match on tax IDs
        nt_match = (czid_s['NT.species_taxid'] == seqtoid_s['NT.species_taxid']).all()
        nr_match = (czid_s['NR.species_taxid'] == seqtoid_s['NR.species_taxid']).all()

        overall_match = nt_match and nr_match

        result_row = {
            'sample': sample,
            'contigs_count_czid': len(czid_s),
            'contigs_count_seqtoid': len(seqtoid_s),
            'NT_species_taxid_match': 'T' if nt_match else 'F',
            'NR_species_taxid_match': 'T' if nr_match else 'F',
            'overall_match': 'T' if overall_match else 'F'
        }

        rows.append(result_row)

    # Save results
    pd.DataFrame(rows).to_csv('step5_contig_summary_reports.csv', index=False)

    if missing_czid or missing_seqtoid:
        pd.DataFrame({
            'missing_czid': [', '.join(sorted(missing_czid))],
            'missing_seqtoid': [', '.join(sorted(missing_seqtoid))]
        }).to_csv('step5_missing.csv', index=False)

    print("  → Results saved to step5_contig_summary_reports.csv")


# Step 6 & 7 remain unchanged (they already handle missing files well)

def compare_nonhost_fastqs():
    print("\n=== Step 6: Non-host reads FASTQ Comparison ===")
    print("    → Checking single unpaired *_reads_nh.fastq per sample")

    rows = []
    missing_czid = []
    missing_seqtoid = []

    for sample in EXPECTED_SAMPLES:
        pattern = f"{sample}_*_reads_nh.fastq"

        czid_files = glob.glob(os.path.join(CZID_DIR, pattern))
        seqtoid_files = glob.glob(os.path.join(SEQTOID_DIR, pattern))

        if len(czid_files) != 1 or len(seqtoid_files) != 1:
            status = 'missing' if len(czid_files) == 0 else 'multiple files'
            if len(czid_files) != 1: missing_czid.append(sample)
            if len(seqtoid_files) != 1: missing_seqtoid.append(sample)
            rows.append({'sample': sample, 'identical': status})
            continue

        czid_fq = czid_files[0]
        seqtoid_fq = seqtoid_files[0]

        # 1. Byte-for-byte identity
        czid_hash = file_sha256(czid_fq)
        seqtoid_hash = file_sha256(seqtoid_fq)

        if czid_hash == seqtoid_hash:
            rows.append({'sample': sample, 'identical': 'T'})
            continue

        row = {'sample': sample, 'identical': 'F'}

        # 2. Line count check (proxy for read count)
        try:
            czid_lines = sum(1 for _ in open(czid_fq))
            seqtoid_lines = sum(1 for _ in open(seqtoid_fq))
            if czid_lines != seqtoid_lines:
                row['line_count'] = f"czid: {czid_lines}, seqtoid: {seqtoid_lines}"
        except Exception as e:
            row['line_count'] = f"error: {str(e)}"

        # 3. Sequence content (sorted, ignore order & qualities)
        try:
            cmd_czid   = f"seqkit seq -s -i {czid_fq}   | sort | sha256sum"
            cmd_seqtoid = f"seqkit seq -s -i {seqtoid_fq} | sort | sha256sum"

            czid_seq_hash   = subprocess.check_output(cmd_czid,   shell=True, text=True).split()[0]
            seqtoid_seq_hash = subprocess.check_output(cmd_seqtoid, shell=True, text=True).split()[0]

            row['sequences_match'] = 'T' if czid_seq_hash == seqtoid_seq_hash else 'F'

        except FileNotFoundError:
            row['sequences_match'] = 'seqkit not found'
        except subprocess.CalledProcessError as e:
            row['sequences_match'] = f"seqkit error: {str(e)}"

        rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv('step6_nonhost_fastq.csv', index=False)

    if missing_czid or missing_seqtoid:
        pd.DataFrame({
            'missing_czid':   [', '.join(sorted(missing_czid))],
            'missing_seqtoid': [', '.join(sorted(missing_seqtoid))]
        }).to_csv('step6_missing.csv', index=False)

    print("  → Results saved to step6_nonhost_fastq.csv")


def compare_nonhost_contigs():
    print("\n=== Step 7: Non-host contigs FASTA Comparison ===")
    rows = []
    missing_czid = []
    missing_seqtoid = []

    for sample in EXPECTED_SAMPLES:
        pattern = f"{sample}_*_contigs_nh.fasta"
        czid_files = glob.glob(os.path.join(CZID_DIR, pattern))
        seqtoid_files = glob.glob(os.path.join(SEQTOID_DIR, pattern))

        if len(czid_files) != 1 or len(seqtoid_files) != 1:
            status = 'missing' if len(czid_files) == 0 else 'multiple files'
            if len(czid_files) != 1: missing_czid.append(sample)
            if len(seqtoid_files) != 1: missing_seqtoid.append(sample)
            rows.append({'sample': sample, 'identical': status})
            continue

        czid_fa = czid_files[0]
        seqtoid_fa = seqtoid_files[0]

        row = {'sample': sample}

        # Primary check: byte-for-byte hash
        czid_hash = file_sha256(czid_fa)
        seqtoid_hash = file_sha256(seqtoid_fa)

        if czid_hash == seqtoid_hash:
            row['identical'] = 'T'
            rows.append(row)
            continue

        row['identical'] = 'F'

        # Secondary check: total sequence length (dependency-free)
        czid_sum = get_fasta_total_length(czid_fa)
        seqtoid_sum = get_fasta_total_length(seqtoid_fa)

        diff = abs(czid_sum - seqtoid_sum)

        row['total_bp_czid']   = czid_sum
        row['total_bp_seqtoid'] = seqtoid_sum
        row['bp_diff']         = diff

        if diff == 0:
            row['length_match'] = 'T (exact)'
        elif diff <= 5:
            row['length_match'] = 'T (≤5 bp)'
        elif diff <= 0.005 * max(czid_sum, seqtoid_sum, 1):
            row['length_match'] = 'equivalent (≤0.005 rel)'
        elif diff <= 0.05 * max(czid_sum, seqtoid_sum, 1):
            row['length_match'] = 'warning (≤0.05 rel)'
        else:
            row['length_match'] = 'significant (>0.05 rel)'

        if czid_sum > 0 and seqtoid_sum > 0:
            rel_diff_pct = (diff / max(czid_sum, seqtoid_sum)) * 100
            row['rel_diff_pct'] = f"{rel_diff_pct:.6f}%"
        else:
            row['rel_diff_pct'] = 'N/A'

        # Optional: sequence content hash (still requires seqkit)
        # Comment out or remove if you want fully dependency-free
        try:
            cmd_czid   = f"seqkit seq -s -i {czid_fa}   | sort | sha256sum"
            cmd_seqtoid = f"seqkit seq -s -i {seqtoid_fa} | sort | sha256sum"

            czid_seq_hash   = subprocess.check_output(cmd_czid,   shell=True, text=True).split()[0]
            seqtoid_seq_hash = subprocess.check_output(cmd_seqtoid, shell=True, text=True).split()[0]

            row['sequences_match'] = 'T' if czid_seq_hash == seqtoid_seq_hash else 'F'

        except Exception as e:
            row['sequences_match'] = f'seqkit unavailable: {str(e)}'

        rows.append(row)

    # Save results
    df = pd.DataFrame(rows)
    df.to_csv('step7_nonhost_contigs.csv', index=False)

    if missing_czid or missing_seqtoid:
        pd.DataFrame({
            'missing_czid':   [', '.join(sorted(set(missing_czid)))],
            'missing_seqtoid': [', '.join(sorted(set(missing_seqtoid)))]
        }).to_csv('step7_missing.csv', index=False)

    print("  → Results saved to step7_nonhost_contigs.csv")


def main():
    print("Starting CZID long-reads pipeline comparison...\n")

    print("=== Step 1: Sample Metadata Comparison ===")
    compare_metadata()

    print("\n=== Step 2: Sample Overviews Comparison ===")
    compare_sample_overviews()

    print("\n=== Step 3: Sample Taxon Reports Comparison ===")
    compare_taxon_reports()

    print("\n=== Step 4: Combined Taxon BPM Comparison ===")
    compare_combined_taxon_results()

    print("\n=== Step 5: Contig Summary Reports Comparison ===")
    compare_contig_summary_reports()

    print("\n=== Step 6: Non-host reads FASTQ Comparison ===")
    compare_nonhost_fastqs()

    print("\n=== Step 7: Non-host contigs FASTA Comparison ===")
    compare_nonhost_contigs()

    print("\nComparison complete. Check CSV files in current directory.")


if __name__ == '__main__':
    main()