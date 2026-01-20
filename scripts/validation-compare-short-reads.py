import pandas as pd
import os
import glob
import numpy as np
from typing import List, Dict

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
    """
    Compare numeric columns of two dataframes (assumed same shape/order after sorting).
    Returns a summary string describing the differences.
    """
    # Identify numeric columns
    num_cols = df1.select_dtypes(include=[np.number]).columns.intersection(df2.columns)

    if len(num_cols) == 0:
        return "No numeric columns to compare"

    results = {}
    max_diffs = {}
    for col in num_cols:
        vals1 = df1[col].fillna(0).values.astype(float)
        vals2 = df2[col].fillna(0).values.astype(float)
        cat = numeric_diff(vals1, vals2, atol=atol)
        results[col] = cat
        max_diffs[col] = np.abs(vals1 - vals2).max()

    # Summarize
    categories = list(results.values())
    if all("identical" in c or "equivalent" in c for c in categories):
        return "Numerically equivalent within tolerance (≤ 0.005)"

    summary_lines = []
    warning_cols = [col for col, cat in results.items() if "warning" in cat]
    sig_cols = [col for col, cat in results.items() if "significant" in cat]

    if warning_cols:
        summary_lines.append(f"Minor differences (0.005–0.05) in {len(warning_cols)} column(s): {', '.join(warning_cols[:3])}")
    if sig_cols:
        summary_lines.append(f"Significant differences (>0.05) in {len(sig_cols)} column(s): {', '.join(sig_cols[:3])}")

    # For wide matrices (Step 4), add cell-level info
    if wide_matrix and 'Taxon Name' in df1.columns:
        # Assume first column is Taxon Name, rest are samples
        sample_cols = [c for c in df1.columns if c != 'Taxon Name']
        cell_diffs = []
        for col in sample_cols:
            vals1 = df1[col].fillna(0).astype(float)
            vals2 = df2[col].fillna(0).astype(float)
            diff = np.abs(vals1 - vals2)
            max_per_sample = diff.max()
            if max_per_sample > 0.05:
                cell_diffs.append((col, max_per_sample, (diff > 0.05).sum()))
        if cell_diffs:
            cell_summary = "\n".join(
                f"  Sample {s}: max diff {d:.4f}, {n} cells >0.05"
                for s, d, n in sorted(cell_diffs, key=lambda x: x[1], reverse=True)[:5]
            )
            summary_lines.append(f"Across samples:\n{cell_summary}")

    if not summary_lines:
        return "Only minor floating-point noise detected"

    return "Numeric differences detected:\n  " + "\n  ".join(summary_lines)


# ────────────────────────────────────────────────────────────────
# Comparison functions
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
        # Could add tolerant check here for numeric columns if desired, but keeping strict for now

    # Missing/extra samples check (unchanged)
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
            if len(czid_files) != 1:
                missing_czid.append(sample)
            if len(seqtoid_files) != 1:
                missing_seqtoid.append(sample)
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
    # (kept strict for now – mostly IDs + some numerics, but can be extended later)
    # Implementation unchanged from previous version
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

    print("\nComparison complete.")


if __name__ == '__main__':
    main()