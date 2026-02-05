import pandas as pd
import os
import glob
import numpy as np
import hashlib
import subprocess
from typing import List, Dict
from scipy.sparse import coo_matrix

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
    'differ':       'F'
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
    """Return dict of numeric column → category."""
    num_cols = df1.select_dtypes(include=[np.number]).columns.intersection(df2.columns)
    results = {}
    for col in num_cols:
        vals1 = df1[col].fillna(0).values.astype(float)
        vals2 = df2[col].fillna(0).values.astype(float)
        cat = numeric_diff(vals1, vals2)
        results[col] = cat
    return results


# ────────────────────────────────────────────────────────────────
# Step functions with CSV output
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

    # Missing / extra
    rows = []
    for name, df in [('czid', czid_df), ('seqtoid', seqtoid_df)]:
        actual = set(df['sample_name'].astype(str))
        missing = ', '.join(sorted(set(EXPECTED_SAMPLES) - actual))
        extra = ', '.join(sorted(actual - set(EXPECTED_SAMPLES)))
        if missing or extra:
            rows.append({'directory': name, 'missing': missing, 'extra': extra})
    if rows:
        pd.DataFrame(rows).to_csv('step1_missing_extra.csv', index=False)


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

    # Missing / extra samples (similar to metadata)
    rows = []
    for name, df in [('czid', czid_df), ('seqtoid', seqtoid_df)]:
        actual = set(df['sample_name'].astype(str))
        missing = ', '.join(sorted(set(EXPECTED_SAMPLES) - actual))
        extra = ', '.join(sorted(actual - set(EXPECTED_SAMPLES)))
        if missing or extra:
            rows.append({'directory': name, 'missing': missing, 'extra': extra})
    if rows:
        pd.DataFrame(rows).to_csv('step2_missing_extra.csv', index=False)


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
            rows.append({'sample': sample, 'strict_identical': 'missing'})
            continue

        czid_path = czid_files[0]
        seqtoid_path = seqtoid_files[0]

        czid_df = pd.read_csv(czid_path, dtype=str)
        seqtoid_df = pd.read_csv(seqtoid_path, dtype=str)

        czid_df['tax_id'] = pd.to_numeric(czid_df['tax_id'], errors='coerce')
        seqtoid_df['tax_id'] = pd.to_numeric(seqtoid_df['tax_id'], errors='coerce')

        czid_s = czid_df.sort_values('tax_id').reset_index(drop=True)
        seqtoid_s = seqtoid_df.sort_values('tax_id').reset_index(drop=True)

        strict_identical = czid_s.equals(seqtoid_s)

        result_row = {'sample': sample, 'strict_identical': 'T' if strict_identical else 'F'}

        if not strict_identical:
            num_results = compare_numeric_dfs(czid_s, seqtoid_s, id_cols=['tax_id'])
            for col, cat in num_results.items():
                result_row[col] = DIFF_SYMBOLS.get(cat, cat)

        rows.append(result_row)

    pd.DataFrame(rows).to_csv('step3_taxon_reports.csv', index=False)

    if missing_czid or missing_seqtoid:
        pd.DataFrame({
            'missing_czid': [', '.join(missing_czid)],
            'missing_seqtoid': [', '.join(missing_seqtoid)]
        }).to_csv('step3_missing.csv', index=False)


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

    # Taxa count and set diff
    taxa_diff = []
    if len(czid_s) != len(seqtoid_s):
        taxa_diff.append({'diff_type': 'row_count', 'czid': len(czid_s), 'seqtoid': len(seqtoid_s)})

    czid_taxa = set(czid_s[sort_col])
    seqtoid_taxa = set(seqtoid_s[sort_col])
    if czid_taxa != seqtoid_taxa:
        extra_czid = len(czid_taxa - seqtoid_taxa)
        extra_seqtoid = len(seqtoid_taxa - czid_taxa)
        taxa_diff.append({'diff_type': 'extra_taxa', 'czid_extra': extra_czid, 'seqtoid_extra': extra_seqtoid})

    if taxa_diff:
        pd.DataFrame(taxa_diff).to_csv('step4_taxa_diff.csv', index=False)


def compare_contig_summary_reports():
    print("\n=== Step 5: Contig Summary Reports Comparison ===")
    rows = []
    missing_czid = []
    missing_seqtoid = []

    for sample in EXPECTED_SAMPLES:
        czid_files = glob.glob(os.path.join(CZID_DIR, f"{sample}_*_contig_summary_report.csv"))
        seqtoid_files = glob.glob(os.path.join(SEQTOID_DIR, f"{sample}_*_contig_summary_report.csv"))

        if len(czid_files) != 1 or len(seqtoid_files) != 1:
            if len(czid_files) != 1: missing_czid.append(sample)
            if len(seqtoid_files) != 1: missing_seqtoid.append(sample)
            rows.append({'sample': sample, 'strict_identical': 'missing'})
            continue

        czid_path = czid_files[0]
        seqtoid_path = seqtoid_files[0]

        czid_df = pd.read_csv(czid_path, dtype=str)
        seqtoid_df = pd.read_csv(seqtoid_path, dtype=str)

        sort_candidates = ['contig_name', 'contig_id', 'contig', 'name', 'id', 'Contig', 'ContigID']
        sort_col = next((c for c in sort_candidates if c in czid_df.columns), None)

        if sort_col:
            czid_df[sort_col] = czid_df[sort_col].astype(str)
            seqtoid_df[sort_col] = seqtoid_df[sort_col].astype(str)
            czid_s = czid_df.sort_values(sort_col).reset_index(drop=True)
            seqtoid_s = seqtoid_df.sort_values(sort_col).reset_index(drop=True)
        else:
            czid_s = czid_df.reset_index(drop=True)
            seqtoid_s = seqtoid_df.reset_index(drop=True)

        strict_identical = czid_s.equals(seqtoid_s)

        result_row = {'sample': sample, 'strict_identical': 'T' if strict_identical else 'F'}

        if not strict_identical:
            num_results = compare_numeric_dfs(czid_s, seqtoid_s, id_cols=[sort_col] if sort_col else [])
            for col, cat in num_results.items():
                result_row[col] = DIFF_SYMBOLS.get(cat, cat)

        rows.append(result_row)

    pd.DataFrame(rows).to_csv('step5_contig_summary_reports.csv', index=False)

    if missing_czid or missing_seqtoid:
        pd.DataFrame({
            'missing_czid': [', '.join(missing_czid)],
            'missing_seqtoid': [', '.join(missing_seqtoid)]
        }).to_csv('step5_missing.csv', index=False)


def compare_nonhost_fastqs():
    print("\n=== Step 6: Non-host reads FASTQ (R1 & R2) ===")
    rows = []
    missing_czid = []
    missing_seqtoid = []

    for sample in EXPECTED_SAMPLES:
        for read in ['R1', 'R2']:
            pattern = f"{sample}_*_reads_nh_{read}.fastq"
            czid_f = glob.glob(os.path.join(CZID_DIR, pattern))
            seqtoid_f = glob.glob(os.path.join(SEQTOID_DIR, pattern))

            if len(czid_f) != 1 or len(seqtoid_f) != 1:
                if len(czid_f) != 1: missing_czid.append(f"{sample} {read}")
                if len(seqtoid_f) != 1: missing_seqtoid.append(f"{sample} {read}")
                rows.append({'sample': sample, 'read': read, 'identical': 'missing'})
                continue

            czid_fa = czid_f[0]
            seqtoid_fa = seqtoid_f[0]

            identical = file_sha256(czid_fa) == file_sha256(seqtoid_fa)
            rows.append({'sample': sample, 'read': read, 'identical': 'T' if identical else 'F'})

    pd.DataFrame(rows).to_csv('step6_nonhost_fastq.csv', index=False)

    if missing_czid or missing_seqtoid:
        pd.DataFrame({
            'missing_czid': [', '.join(missing_czid)],
            'missing_seqtoid': [', '.join(missing_seqtoid)]
        }).to_csv('step6_missing.csv', index=False)


def compare_nonhost_contigs():
    print("\n=== Step 7: Non-host contigs FASTA ===")
    rows = []
    missing_czid = []
    missing_seqtoid = []

    for sample in EXPECTED_SAMPLES:
        pattern = f"{sample}_*_contigs_nh.fasta"
        czid_f = glob.glob(os.path.join(CZID_DIR, pattern))
        seqtoid_f = glob.glob(os.path.join(SEQTOID_DIR, pattern))

        if len(czid_f) != 1 or len(seqtoid_f) != 1:
            if len(czid_f) != 1: missing_czid.append(sample)
            if len(seqtoid_f) != 1: missing_seqtoid.append(sample)
            rows.append({'sample': sample, 'identical': 'missing'})
            continue

        czid_fa = czid_f[0]
        seqtoid_fa = seqtoid_f[0]

        identical = file_sha256(czid_fa) == file_sha256(seqtoid_fa)
        rows.append({'sample': sample, 'identical': 'T' if identical else 'F'})

    pd.DataFrame(rows).to_csv('step7_nonhost_contigs.csv', index=False)

    if missing_czid or missing_seqtoid:
        pd.DataFrame({
            'missing_in_czid': [', '.join(missing_czid)],
            'missing_in_seqtoid': [', '.join(missing_seqtoid)]
        }).to_csv('step7_missing.csv', index=False)


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