import pandas as pd
import os
import glob
import hashlib
import subprocess
from typing import List

# ────────────────────────────────────────────────────────────────
# Configuration
# ────────────────────────────────────────────────────────────────

CZID_DIR = 'czid'
SEQTOID_DIR = 'seqtoid'

METADATA_FILE = 'sample_metadata.csv'

# Load expected samples directly from the provided metadata
try:
    metadata_df = pd.read_csv(METADATA_FILE)
    EXPECTED_SAMPLES: List[str] = sorted(metadata_df['sample_name'].astype(str).unique().tolist())
except FileNotFoundError:
    print(f"Error: {METADATA_FILE} not found. Using fallback empty list.")
    EXPECTED_SAMPLES: List[str] = []

# If you want to hard-code for safety/debugging (uncomment if needed):
# EXPECTED_SAMPLES = [
#     'ERR11417004',
#     'SRR10903401',
#     'SRR12876565',
#     'SRR13227003',
#     'SRR13227004',
#     'SRR13227005',
#     'SRR15049352'
# ]

CONTIG_PATTERN = '*_contigs.fa'          # e.g. ERR11417004_881252_contigs.fa

DIFF_SYMBOLS = {
    'equivalent':   '✅ <0.005',
    'warning':      '⚠️ 0.005–0.05',
    'significant':  '❌ >0.05',
    'identical':    'T',
    'differ':       'F',
    'missing':      'missing',
    'content_identical': 'sequences identical',
    'content_diff': 'sequences differ',
}

# ────────────────────────────────────────────────────────────────
# Helpers
# ────────────────────────────────────────────────────────────────

def file_sha256(filepath: str) -> str:
    """Compute SHA-256 hash of a file in chunks."""
    sha256 = hashlib.sha256()
    with open(filepath, 'rb') as f:
        while chunk := f.read(8192 * 1024):
            sha256.update(chunk)
    return sha256.hexdigest()


def compare_numeric_dfs(df1: pd.DataFrame, df2: pd.DataFrame, id_cols: List[str]) -> dict:
    num_cols = df1.select_dtypes(include=['number']).columns.intersection(df2.columns)
    results = {}
    for col in num_cols:
        vals1 = df1[col].fillna(0).values.astype(float)
        vals2 = df2[col].fillna(0).values.astype(float)
        diff = abs(vals1 - vals2).max()
        if pd.isna(diff):
            results[col] = 'NaN'
        elif diff <= 0.005:
            results[col] = 'equivalent'
        elif diff <= 0.05:
            results[col] = 'warning'
        else:
            results[col] = 'significant'
    return results


# ────────────────────────────────────────────────────────────────
# Step 1: Metadata
# ────────────────────────────────────────────────────────────────

def compare_metadata():
    czid_path = os.path.join(CZID_DIR, METADATA_FILE)
    seqtoid_path = os.path.join(SEQTOID_DIR, METADATA_FILE)

    if not os.path.exists(czid_path) or not os.path.exists(seqtoid_path):
        print(f"Error: {METADATA_FILE} missing in one or both directories.")
        pd.DataFrame({'file': [METADATA_FILE], 'identical': ['missing']}).to_csv('step1_metadata.csv', index=False)
        return

    czid_df = pd.read_csv(czid_path)
    seqtoid_df = pd.read_csv(seqtoid_path)

    czid_s = czid_df.sort_values('sample_name').reset_index(drop=True)
    seqtoid_s = seqtoid_df.sort_values('sample_name').reset_index(drop=True)

    identical = czid_s.equals(seqtoid_s)
    pd.DataFrame({'file': [METADATA_FILE], 'identical': ['T' if identical else 'F']}).to_csv('step1_metadata.csv', index=False)

    # Check for missing/extra samples
    rows = []
    for name, df in [('czid', czid_df), ('seqtoid', seqtoid_df)]:
        actual = set(df['sample_name'].astype(str))
        missing = ', '.join(sorted(set(EXPECTED_SAMPLES) - actual))
        extra = ', '.join(sorted(actual - set(EXPECTED_SAMPLES)))
        if missing or extra:
            rows.append({'directory': name, 'missing': missing, 'extra': extra})
    if rows:
        pd.DataFrame(rows).to_csv('step1_missing_extra.csv', index=False)


# ────────────────────────────────────────────────────────────────
# Step 2: Combined AMR Results (gene presence comparison)
# ────────────────────────────────────────────────────────────────

def compare_combined_amr_results():
    file_name = 'combined_amr_results.csv'
    czid_path = os.path.join(CZID_DIR, file_name)
    seqtoid_path = os.path.join(SEQTOID_DIR, file_name)

    print("Step 2: Comparing combined_amr_results.csv (gene presence)...")

    if not os.path.exists(czid_path) or not os.path.exists(seqtoid_path):
        print(f"Error: {file_name} missing in one or both directories.")
        pd.DataFrame({'file': [file_name], 'status': ['missing']}).to_csv('step2_amr_comparison.csv', index=False)
        return

    czid_df = pd.read_csv(czid_path)
    seqtoid_df = pd.read_csv(seqtoid_path)

    results = []
    missing_samples_czid = []
    missing_samples_seqtoid = []

    for sample in EXPECTED_SAMPLES:
        czid_sample = czid_df[czid_df['sample_name'] == sample]
        seqtoid_sample = seqtoid_df[seqtoid_df['sample_name'] == sample]

        if czid_sample.empty:
            missing_samples_czid.append(sample)
        if seqtoid_sample.empty:
            missing_samples_seqtoid.append(sample)

        if czid_sample.empty or seqtoid_sample.empty:
            results.append({
                'sample': sample,
                'missing_proportion': 'N/A',
                'category': 'missing',
                'symbol': 'missing'
            })
            continue

        czid_genes = set(czid_sample['gene_name'].dropna().unique())
        seqtoid_genes = set(seqtoid_sample['gene_name'].dropna().unique())

        missing_in_seqtoid_count = len(czid_genes - seqtoid_genes)
        total_czid_genes = len(czid_genes)

        proportion_missing = missing_in_seqtoid_count / total_czid_genes if total_czid_genes > 0 else 0.0

        if proportion_missing < 0.005:
            category = 'equivalent'
        elif proportion_missing < 0.05:
            category = 'warning'
        else:
            category = 'significant'

        symbol = DIFF_SYMBOLS.get(category, category)

        results.append({
            'sample': sample,
            'total_czid_genes': total_czid_genes,
            'missing_in_seqtoid': missing_in_seqtoid_count,
            'missing_proportion': round(proportion_missing, 6),
            'category': category,
            'symbol': symbol
        })

        print(f"  - {sample}: missing proportion = {proportion_missing:.4f} ({missing_in_seqtoid_count}/{total_czid_genes}) → {symbol}")

    pd.DataFrame(results).to_csv('step2_amr_comparison.csv', index=False)

    if missing_samples_czid or missing_samples_seqtoid:
        pd.DataFrame({
            'missing_in_czid': [', '.join(sorted(missing_samples_czid))],
            'missing_in_seqtoid': [', '.join(sorted(missing_samples_seqtoid))]
        }).to_csv('step2_amr_missing.csv', index=False)


# ────────────────────────────────────────────────────────────────
# Step 3: Contig FASTA Comparison
# ────────────────────────────────────────────────────────────────

def compare_contig_fastas():
    print("Step 3: Comparing non-host contigs FASTA files (*_contigs.fa)...")

    rows = []
    missing_czid = []
    missing_seqtoid = []

    for sample in EXPECTED_SAMPLES:
        czid_files = glob.glob(os.path.join(CZID_DIR, f"{sample}{CONTIG_PATTERN}"))
        seqtoid_files = glob.glob(os.path.join(SEQTOID_DIR, f"{sample}{CONTIG_PATTERN}"))

        if len(czid_files) != 1 or len(seqtoid_files) != 1:
            if len(czid_files) != 1:
                missing_czid.append(sample)
            if len(seqtoid_files) != 1:
                missing_seqtoid.append(sample)
            rows.append({'sample': sample, 'identical': 'missing', 'note': ''})
            print(f"  - {sample}: missing")
            continue

        czid_path = czid_files[0]
        seqtoid_path = seqtoid_files[0]

        file_hash_identical = file_sha256(czid_path) == file_sha256(seqtoid_path)

        if file_hash_identical:
            status = 'identical'
            note = ''
        else:
            # Deeper check: compare sorted sequences only (ignores order & headers)
            try:
                cmd = (
                    f"seqkit seq -s '{czid_path}' | sort | sha256sum | cut -d' ' -f1 && "
                    f"seqkit seq -s '{seqtoid_path}' | sort | sha256sum | cut -d' ' -f1"
                )
                result = subprocess.check_output(cmd, shell=True, text=True).strip().split('\n')
                if len(result) == 2 and result[0] == result[1]:
                    status = 'content_identical'
                    note = '(sequences identical, headers/order differ)'
                else:
                    status = 'content_diff'
                    note = '(sequences differ)'
            except Exception:
                status = 'differ'
                note = '(content check failed)'

        rows.append({
            'sample': sample,
            'identical': DIFF_SYMBOLS.get(status, status),
            'note': note
        })

        print(f"  - {sample}: {DIFF_SYMBOLS.get(status, status)} {note}")

    pd.DataFrame(rows).to_csv('step3_contigs_fasta.csv', index=False)

    if missing_czid or missing_seqtoid:
        pd.DataFrame({
            'missing_in_czid': [', '.join(sorted(missing_czid))],
            'missing_in_seqtoid': [', '.join(sorted(missing_seqtoid))]
        }).to_csv('step3_missing.csv', index=False)


# ────────────────────────────────────────────────────────────────
# Main
# ────────────────────────────────────────────────────────────────

def main():
    print("Starting CZID AMR pipeline comparison...\n")
    print(f"Expected samples ({len(EXPECTED_SAMPLES)}): {', '.join(EXPECTED_SAMPLES)}\n")

    print("=== Step 1: Sample Metadata Comparison ===")
    compare_metadata()

    print("\n=== Step 2: Combined AMR Results (gene presence) Comparison ===")
    compare_combined_amr_results()

    print("\n=== Step 3: Non-host Contigs FASTA Comparison ===")
    compare_contig_fastas()

    print("\nComparison complete. Check CSV files in current directory.")


if __name__ == '__main__':
    main()