import pandas as pd
import os
import glob
import numpy as np
import json
from typing import List, Dict
from scipy.sparse import coo_matrix

# Define paths to the directories
CZID_DIR = 'czid'
SEQTOID_DIR = 'seqtoid'

# List of expected sample IDs
EXPECTED_SAMPLES = [
    'ERR11417004',
    'SRR23038836',
    'SRR11454628',
    'SRR12876565',
    'SRR14579537',
    'SRR10903401',
    'SRR1304850',
    'SRR13227005',
    'SRR13227004',
    'SRR13227003',
    'SRR11278904'
]

# ────────────────────────────────────────────────────────────────
# Helper functions for tolerant numeric comparison
# ────────────────────────────────────────────────────────────────

def numeric_diff(a: np.ndarray, b: np.ndarray, atol: float = 0.005) -> str:
    """Categorize the worst difference between two numeric arrays."""
    if len(a) == 0 or len(b) == 0:
        return "empty arrays"
    diff = np.abs(a - b)
    worst = diff.max()
    if np.isnan(worst):
        return "contains NaN"
    if worst <= 1e-8:
        return "identical"
    elif worst <= 0.005:
        return f"equivalent (max diff {worst:.6f} ≤ 0.005)"
    elif worst <= 0.05:
        return f"warning (max diff {worst:.6f})"
    else:
        return f"significant (max diff {worst:.6f})"


def compare_numeric_dfs(
        df1: pd.DataFrame,
        df2: pd.DataFrame,
        id_cols: List[str],
        atol: float = 0.005,
        wide_matrix: bool = False
) -> str:
    """Compare numeric columns of two dataframes with tolerance."""
    num_cols = df1.select_dtypes(include=[np.number]).columns.intersection(df2.columns)
    if len(num_cols) == 0:
        return "No numeric columns to compare"

    results = {}
    for col in num_cols:
        vals1 = df1[col].fillna(0).values.astype(float)
        vals2 = df2[col].fillna(0).values.astype(float)
        cat = numeric_diff(vals1, vals2, atol=atol)
        results[col] = cat

    categories = list(results.values())
    if all("identical" in c or "equivalent" in c for c in categories):
        return "Numerically equivalent within tolerance (≤ 0.005)"

    summary_lines = []
    warning_cols = [col for col, cat in results.items() if "warning" in cat]
    sig_cols = [col for col, cat in results.items() if "significant" in cat]

    if warning_cols:
        summary_lines.append(f"Minor differences (0.005–0.05) in {len(warning_cols)} column(s)")
    if sig_cols:
        summary_lines.append(f"Significant differences (>0.05) in {len(sig_cols)} column(s)")

    if not summary_lines:
        return "Only minor floating-point noise detected"

    return "Numeric differences detected:\n  " + "\n  ".join(summary_lines)


# ────────────────────────────────────────────────────────────────
# Comparison functions (Steps 1–6)
# ────────────────────────────────────────────────────────────────

def compare_metadata():
    metadata_file = 'sample_metadata.csv'
    czid_path = os.path.join(CZID_DIR, metadata_file)
    seqtoid_path = os.path.join(SEQTOID_DIR, metadata_file)

    if not os.path.exists(czid_path) or not os.path.exists(seqtoid_path):
        print(f"Error: One or both {metadata_file} files missing.")
        return

    czid_df = pd.read_csv(czid_path)
    seqtoid_df = pd.read_csv(seqtoid_path)

    czid_sorted = czid_df.sort_values('sample_name').reset_index(drop=True)
    seqtoid_sorted = seqtoid_df.sort_values('sample_name').reset_index(drop=True)

    print("Comparing sample_metadata.csv...")
    if czid_sorted.equals(seqtoid_sorted):
        print("  - Files are identical.")
    else:
        print("  - Files differ (strict comparison failed).")

    for name, df in [('czid', czid_df), ('seqtoid', seqtoid_df)]:
        actual = set(df['sample_name'].astype(str))
        expected = set(EXPECTED_SAMPLES)
        missing = expected - actual
        extra = actual - expected
        if missing:
            print(f"  - Missing in {name}: {', '.join(sorted(missing))}")
        if extra:
            print(f"  - Extra in {name}: {', '.join(sorted(extra))}")


def compare_overviews():
    file_name = 'sample_overviews.csv'
    czid_path = os.path.join(CZID_DIR, file_name)
    seqtoid_path = os.path.join(SEQTOID_DIR, file_name)

    if not os.path.exists(czid_path) or not os.path.exists(seqtoid_path):
        print(f"Error: One or both {file_name} files missing.")
        return

    czid_df = pd.read_csv(czid_path, dtype=str)
    seqtoid_df = pd.read_csv(seqtoid_path, dtype=str)

    czid_sorted = czid_df.sort_values('sample_name').reset_index(drop=True)
    seqtoid_sorted = seqtoid_df.sort_values('sample_name').reset_index(drop=True)

    print("Comparing sample_overviews.csv...")
    if czid_sorted.equals(seqtoid_sorted):
        print("  - Files are identical.")
    else:
        print("  - Strict comparison failed → checking numeric tolerance...")
        result = compare_numeric_dfs(czid_sorted, seqtoid_sorted, id_cols=['sample_name'])
        print("    → " + result)


def compare_taxon_reports():
    print("Comparing per-sample taxon_report.csv files...")
    missing_czid = []
    missing_seqtoid = []

    for sample in EXPECTED_SAMPLES:
        czid_files = glob.glob(os.path.join(CZID_DIR, f"{sample}_*_taxon_report.csv"))
        seqtoid_files = glob.glob(os.path.join(SEQTOID_DIR, f"{sample}_*_taxon_report.csv"))

        if len(czid_files) != 1 or len(seqtoid_files) != 1:
            if len(czid_files) != 1: missing_czid.append(sample)
            if len(seqtoid_files) != 1: missing_seqtoid.append(sample)
            continue

        czid_df = pd.read_csv(czid_files[0], dtype=str)
        seqtoid_df = pd.read_csv(seqtoid_files[0], dtype=str)

        czid_df['tax_id'] = pd.to_numeric(czid_df['tax_id'], errors='coerce')
        seqtoid_df['tax_id'] = pd.to_numeric(seqtoid_df['tax_id'], errors='coerce')

        czid_sorted = czid_df.sort_values('tax_id').reset_index(drop=True)
        seqtoid_sorted = seqtoid_df.sort_values('tax_id').reset_index(drop=True)

        print(f"  - {sample} ...", end=" ")
        if czid_sorted.equals(seqtoid_sorted):
            print("identical")
        else:
            print("strict fail → tolerant check")
            result = compare_numeric_dfs(czid_sorted, seqtoid_sorted, id_cols=['tax_id'])
            print(f"    → {result}")

    if missing_czid:
        print(f"  Missing in czid: {', '.join(missing_czid)}")
    if missing_seqtoid:
        print(f"  Missing in seqtoid: {', '.join(missing_seqtoid)}")


def compare_combined_taxon_results():
    file_name = 'combined_sample_taxon_results_NT.rpm.csv'
    czid_path = os.path.join(CZID_DIR, file_name)
    seqtoid_path = os.path.join(SEQTOID_DIR, file_name)

    if not os.path.exists(czid_path) or not os.path.exists(seqtoid_path):
        print(f"Error: One or both {file_name} files missing.")
        return

    czid_df = pd.read_csv(czid_path, dtype=str)
    seqtoid_df = pd.read_csv(seqtoid_path, dtype=str)

    sort_col = 'Taxon Name' if 'Taxon Name' in czid_df.columns else czid_df.columns[0]
    czid_sorted = czid_df.sort_values(sort_col).reset_index(drop=True)
    seqtoid_sorted = seqtoid_df.sort_values(sort_col).reset_index(drop=True)

    print("Comparing combined_sample_taxon_results_NT.rpm.csv...")
    if czid_sorted.equals(seqtoid_sorted):
        print("  - Files are identical.")
    else:
        print("  - Strict comparison failed → checking numeric tolerance (wide matrix)...")
        result = compare_numeric_dfs(czid_sorted, seqtoid_sorted, id_cols=[sort_col], wide_matrix=True)
        print("    → " + result)


def compare_contig_summary_reports():
    print("Comparing per-sample contig_summary_report.csv files...")
    missing_czid = []
    missing_seqtoid = []

    for sample in EXPECTED_SAMPLES:
        czid_files = glob.glob(os.path.join(CZID_DIR, f"{sample}_*_contig_summary_report.csv"))
        seqtoid_files = glob.glob(os.path.join(SEQTOID_DIR, f"{sample}_*_contig_summary_report.csv"))

        if len(czid_files) != 1 or len(seqtoid_files) != 1:
            if len(czid_files) != 1: missing_czid.append(sample)
            if len(seqtoid_files) != 1: missing_seqtoid.append(sample)
            continue

        czid_df = pd.read_csv(czid_files[0], dtype=str)
        seqtoid_df = pd.read_csv(seqtoid_files[0], dtype=str)

        czid_sorted = czid_df.sort_values('contig_name').reset_index(drop=True)
        seqtoid_sorted = seqtoid_df.sort_values('contig_name').reset_index(drop=True)

        print(f"  - {sample} ...", end=" ")
        if czid_sorted.equals(seqtoid_sorted):
            print("identical")
        else:
            print("differ (strict comparison)")

    if missing_czid:
        print(f"  Missing in czid: {', '.join(missing_czid)}")
    if missing_seqtoid:
        print(f"  Missing in seqtoid: {', '.join(missing_seqtoid)}")


def compare_host_gene_counts():
    print("Comparing per-sample reads_per_transcript.kallisto.tsv files...")
    missing_czid = []
    missing_seqtoid = []

    for sample in EXPECTED_SAMPLES:
        czid_files = glob.glob(os.path.join(CZID_DIR, f"{sample}_*_reads_per_transcript.kallisto.tsv"))
        seqtoid_files = glob.glob(os.path.join(SEQTOID_DIR, f"{sample}_*_reads_per_transcript.kallisto.tsv"))

        if len(czid_files) != 1 or len(seqtoid_files) != 1:
            if len(czid_files) != 1: missing_czid.append(sample)
            if len(seqtoid_files) != 1: missing_seqtoid.append(sample)
            continue

        czid_df = pd.read_csv(czid_files[0], sep='\t', dtype=str)
        seqtoid_df = pd.read_csv(seqtoid_files[0], sep='\t', dtype=str)

        czid_sorted = czid_df.sort_values('target_id').reset_index(drop=True)
        seqtoid_sorted = seqtoid_df.sort_values('target_id').reset_index(drop=True)

        print(f"  - {sample} ...", end=" ")
        if czid_sorted.equals(seqtoid_sorted):
            print("identical")
        else:
            print("strict fail → tolerant check")
            result = compare_numeric_dfs(czid_sorted, seqtoid_sorted, id_cols=['target_id'])
            print(f"    → {result}")

    if missing_czid:
        print(f"  Missing in czid: {', '.join(missing_czid)}")
    if missing_seqtoid:
        print(f"  Missing in seqtoid: {', '.join(missing_seqtoid)}")


# ────────────────────────────────────────────────────────────────
# Step 7: Combined Microbiome File (BIOM) with sparse matrix diff
# ────────────────────────────────────────────────────────────────

def compare_combined_microbiome():
    biom_file = 'Combined Microbiome File.biom'  # adjust filename if different

    czid_path = os.path.join(CZID_DIR, biom_file)
    seqtoid_path = os.path.join(SEQTOID_DIR, biom_file)

    if not os.path.exists(czid_path):
        print(f"Error: {czid_path} not found.")
        return
    if not os.path.exists(seqtoid_path):
        print(f"Error: {seqtoid_path} not found.")
        return

    with open(czid_path, 'r') as f:
        czid_biom = json.load(f)
    with open(seqtoid_path, 'r') as f:
        seqtoid_biom = json.load(f)

    print("Comparing Combined Microbiome File.biom...")

    # Basic structure check
    keys_match = czid_biom.keys() == seqtoid_biom.keys()
    shape_match = czid_biom.get('shape') == seqtoid_biom.get('shape')
    matrix_type_match = czid_biom.get('matrix_type') == seqtoid_biom.get('matrix_type')

    if not (keys_match and shape_match and matrix_type_match):
        print("  - Basic structure differs:")
        if not keys_match: print("    Different top-level keys")
        if not shape_match: print(f"    Shape: czid {czid_biom.get('shape')}, seqtoid {seqtoid_biom.get('shape')}")
        if not matrix_type_match: print("    Matrix type differs")
        return

    print(f"  - Shape: {czid_biom['shape']} (rows=features/taxa, cols=samples)")

    # Rows (taxa)
    czid_rows = sorted(czid_biom['rows'], key=lambda x: x['id'])
    seqtoid_rows = sorted(seqtoid_biom['rows'], key=lambda x: x['id'])
    print(f"  - Taxa/rows: {'identical' if czid_rows == seqtoid_rows else 'differ'}")

    # Columns (samples)
    czid_cols = sorted(czid_biom['columns'], key=lambda x: x['id'])
    seqtoid_cols = sorted(seqtoid_biom['columns'], key=lambda x: x['id'])
    print(f"  - Samples/columns: {'identical' if czid_cols == seqtoid_cols else 'differ'}")

    # Sample name check (strip serial)
    czid_sample_ids = {col['id'].split(':')[0] for col in czid_biom['columns']}
    seqtoid_sample_ids = {col['id'].split(':')[0] for col in seqtoid_biom['columns']}
    expected_set = set(EXPECTED_SAMPLES)
    for name, actual in [('czid', czid_sample_ids), ('seqtoid', seqtoid_sample_ids)]:
        missing = expected_set - actual
        extra = actual - expected_set
        if missing:
            print(f"  - Missing samples in {name}: {', '.join(sorted(missing))}")
        if extra:
            print(f"  - Extra samples in {name}: {', '.join(sorted(extra))}")

    # Sparse data with coo_matrix
    czid_data = sorted(czid_biom['data'], key=lambda x: (x[0], x[1]))
    seqtoid_data = sorted(seqtoid_biom['data'], key=lambda x: (x[0], x[1]))

    if len(czid_data) != len(seqtoid_data):
        print(f"  - Non-zero count differs: czid {len(czid_data)}, seqtoid {len(seqtoid_data)}")
        return

    czid_pos = [(x[0], x[1]) for x in czid_data]
    seqtoid_pos = [(x[0], x[1]) for x in seqtoid_data]

    if czid_pos != seqtoid_pos:
        print("  - Sparsity pattern differs (different positions have non-zero values)")
        return

    czid_values = np.array([x[2] for x in czid_data], dtype=float)
    seqtoid_values = np.array([x[2] for x in seqtoid_data], dtype=float)

    print("  - Sparse positions match → comparing abundance values...")
    cat = numeric_diff(czid_values, seqtoid_values)
    print(f"    → {cat}")

    # Build sparse matrices for extra validation
    shape = tuple(czid_biom['shape'])
    czid_sparse = coo_matrix(
        (czid_values, (np.array([x[0] for x in czid_data]), np.array([x[1] for x in czid_data]))),
        shape=shape
    )
    seqtoid_sparse = coo_matrix(
        (seqtoid_values, (np.array([x[0] for x in seqtoid_data]), np.array([x[1] for x in seqtoid_data]))),
        shape=shape
    )

    max_abs_diff = np.abs(czid_sparse.data - seqtoid_sparse.data).max()
    mean_abs_diff = np.abs(czid_sparse.data - seqtoid_sparse.data).mean()
    nz_count = len(czid_sparse.data)

    print(f"    Non-zero elements: {nz_count}")
    print(f"    Max absolute difference: {max_abs_diff:.6f}")
    print(f"    Mean absolute difference: {mean_abs_diff:.6f}")


def main():
    print("Starting CZID pipeline output comparison...\n")

    print("=== Step 1: Sample Metadata Comparison ===")
    compare_metadata()

    print("\n=== Step 2: Sample Overviews Comparison ===")
    compare_overviews()

    print("\n=== Step 3: Sample Taxon Reports Comparison ===")
    compare_taxon_reports()

    print("\n=== Step 4: Combined Sample Taxon Results Comparison ===")
    compare_combined_taxon_results()

    print("\n=== Step 5: Sample Contig Summary Reports Comparison ===")
    compare_contig_summary_reports()

    print("\n=== Step 6: Sample Host Gene Counts Comparison ===")
    compare_host_gene_counts()

    print("\n=== Step 7: Combined Microbiome File (BIOM) Comparison ===")
    compare_combined_microbiome()

    print("\nComparison complete.")


if __name__ == '__main__':
    main()