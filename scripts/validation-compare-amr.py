import pandas as pd
import os
from typing import List

# ────────────────────────────────────────────────────────────────
# Configuration
# ────────────────────────────────────────────────────────────────

CZID_DIR = 'czid'
SEQTOID_DIR = 'seqtoid'

# Load expected samples from metadata
METADATA_FILE = 'sample_metadata.csv'
metadata_df = pd.read_csv(METADATA_FILE)
EXPECTED_SAMPLES: List[str] = metadata_df['sample_name'].tolist()

DIFF_SYMBOLS = {
    'equivalent':   '✅ <0.005',
    'warning':      '⚠️ 0.005–0.05',
    'significant':  '❌ >0.05',
    'identical':    'T',
    'differ':       'F'
}

# ────────────────────────────────────────────────────────────────
# Helpers
# ────────────────────────────────────────────────────────────────

def compare_metadata():
    czid_path = os.path.join(CZID_DIR, METADATA_FILE)
    seqtoid_path = os.path.join(SEQTOID_DIR, METADATA_FILE)

    if not os.path.exists(czid_path) or not os.path.exists(seqtoid_path):
        print(f"Error: {METADATA_FILE} missing in one or both dirs.")
        pd.DataFrame({'file': [METADATA_FILE], 'identical': ['missing']}).to_csv('step1_metadata.csv', index=False)
        return

    czid_df = pd.read_csv(czid_path)
    seqtoid_df = pd.read_csv(seqtoid_path)

    czid_s = czid_df.sort_values('sample_name').reset_index(drop=True)
    seqtoid_s = seqtoid_df.sort_values('sample_name').reset_index(drop=True)

    identical = czid_s.equals(seqtoid_s)
    pd.DataFrame({'file': [METADATA_FILE], 'identical': ['T' if identical else 'F']}).to_csv('step1_metadata.csv', index=False)

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

# ────────────────────────────────────────────────────────────────
# Step 2: Combined AMR Results Comparison
# ────────────────────────────────────────────────────────────────

def compare_combined_amr_results():
    """
    Step 2: Compare combined_amr_results.csv files.
    - For each sample, genes in czid but not in seqtoid.
    - Proportion missing = (czid_genes - seqtoid_genes) / total_czid_genes
    - Categories: <0.005 equivalent, 0.005-0.05 warning, >0.05 significant.
    - Writes per-sample results to CSV.
    """
    file_name = 'combined_amr_results.csv'
    czid_path = os.path.join(CZID_DIR, file_name)
    seqtoid_path = os.path.join(SEQTOID_DIR, file_name)

    print("Comparing combined_amr_results.csv...")

    if not os.path.exists(czid_path) or not os.path.exists(seqtoid_path):
        print(f"Error: One or both {file_name} files missing.")
        pd.DataFrame({'file': [file_name], 'status': ['missing']}).to_csv('step2_amr_comparison.csv', index=False)
        return

    czid_df = pd.read_csv(czid_path)
    seqtoid_df = pd.read_csv(seqtoid_path)

    results = []
    missing_in_czid = []
    missing_in_seqtoid = []

    for sample in EXPECTED_SAMPLES:
        # Filter for sample
        czid_sample_df = czid_df[czid_df['sample_name'] == sample]
        seqtoid_sample_df = seqtoid_df[seqtoid_df['sample_name'] == sample]

        if czid_sample_df.empty:
            missing_in_czid.append(sample)
        if seqtoid_sample_df.empty:
            missing_in_seqtoid.append(sample)

        if czid_sample_df.empty or seqtoid_sample_df.empty:
            results.append({'sample': sample, 'missing_proportion': 'N/A', 'category': 'missing', 'symbol': 'missing'})
            continue

        # Get unique gene sets (using gene_name)
        czid_genes = set(czid_sample_df['gene_name'].unique())
        seqtoid_genes = set(seqtoid_sample_df['gene_name'].unique())

        # Counts
        missing_in_seqtoid = len(czid_genes - seqtoid_genes)
        total_czid = len(czid_genes)

        # Proportion of czid genes missing in seqtoid
        missing_proportion = missing_in_seqtoid / total_czid if total_czid > 0 else 0

        # Categorize
        if missing_proportion < 0.005:
            category = 'equivalent'
        elif missing_proportion < 0.05:
            category = 'warning'
        else:
            category = 'significant'

        symbol = DIFF_SYMBOLS.get(category, category)

        print(f"  - {sample}: missing proportion {missing_proportion:.4f} (missing {missing_in_seqtoid}, total czid {total_czid}) → {symbol}")

        results.append({
            'sample': sample,
            'missing_in_seqtoid': missing_in_seqtoid,
            'total_czid_genes': total_czid,
            'missing_proportion': round(missing_proportion, 6),
            'category': category,
            'symbol': symbol
        })

    # Write CSV
    pd.DataFrame(results).to_csv('step2_amr_comparison.csv', index=False)

    if missing_in_czid or missing_in_seqtoid:
        pd.DataFrame({
            'missing_in_czid': [', '.join(missing_in_czid)],
            'missing_in_seqtoid': [', '.join(missing_in_seqtoid)]
        }).to_csv('step2_missing.csv', index=False)

# ────────────────────────────────────────────────────────────────
# Main
# ────────────────────────────────────────────────────────────────

def main():
    print("Starting CZID AMR pipeline comparison...\n")

    print("=== Step 1: Sample Metadata Comparison ===")
    compare_metadata()

    print("\n=== Step 2: Combined AMR Results Comparison ===")
    compare_combined_amr_results()

    print("\nComparison complete. Check CSV files in current directory.")

if __name__ == '__main__':
    main()