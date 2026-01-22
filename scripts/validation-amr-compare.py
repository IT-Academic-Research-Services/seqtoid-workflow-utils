#!/usr/bin/env python3
"""
CZID AMR Pipeline Comparison – Focused on core downloadable files
Compares: sample_metadata.csv + combined_amr_results.csv

Run from directory containing:
    czid/
    seqtoid/
"""

import os
import pandas as pd
import numpy as np
import hashlib
from typing import List

# ────────────────────────────────────────────────────────────────
# Configuration
# ────────────────────────────────────────────────────────────────

CZID_DIR    = 'czid'
SEQTOID_DIR = 'seqtoid'

EXPECTED_SAMPLES: List[str] = [
    'ERR11417004',
    'SRR10903401',
    'SRR12876565',
    'SRR13227003',
    'SRR13227004',
    'SRR13227005',
    'SRR15049352',
    # 'SRR14579537_75M'  # absent from metadata & combined file
]

# ────────────────────────────────────────────────────────────────
# Shared helpers
# ────────────────────────────────────────────────────────────────

def file_sha256(filepath: str) -> str:
    sha256 = hashlib.sha256()
    with open(filepath, 'rb') as f:
        while chunk := f.read(8192 * 1024):
            sha256.update(chunk)
    return sha256.hexdigest()


def numeric_diff(a: np.ndarray, b: np.ndarray, atol: float = 0.005) -> str:
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
        id_cols: list = None,
        atol: float = 0.005,
) -> str:
    if id_cols is None:
        id_cols = []
    num_cols = df1.select_dtypes(include=[np.number]).columns.intersection(df2.columns)
    if len(num_cols) == 0:
        return "No numeric columns found to compare"

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
    sig_cols     = [col for col, cat in results.items() if "significant" in cat]

    if warning_cols:
        summary_lines.append(f"Minor differences (0.005–0.05) in {len(warning_cols)} column(s)")
    if sig_cols:
        summary_lines.append(f"Significant differences (>0.05) in {len(sig_cols)} column(s)")

    if not summary_lines:
        return "Only minor floating-point noise detected"

    return "Numeric differences detected:\n  " + "\n  ".join(summary_lines)


# ────────────────────────────────────────────────────────────────
# Step 1: Sample Metadata
# ────────────────────────────────────────────────────────────────

def compare_metadata():
    print("\n=== Step 1: Sample Metadata Comparison ===")
    metadata_file = 'sample_metadata.csv'
    czid_path = os.path.join(CZID_DIR, metadata_file)
    seqtoid_path = os.path.join(SEQTOID_DIR, metadata_file)

    if not all(os.path.exists(p) for p in [czid_path, seqtoid_path]):
        print("✗ One or both sample_metadata.csv files missing.")
        return

    czid_df    = pd.read_csv(czid_path)
    seqtoid_df = pd.read_csv(seqtoid_path)

    czid_sorted    = czid_df.sort_values('sample_name').reset_index(drop=True)
    seqtoid_sorted = seqtoid_df.sort_values('sample_name').reset_index(drop=True)

    if czid_sorted.equals(seqtoid_sorted):
        print("  ✓ Files are identical (strict equality after sorting)")
    else:
        print("  ⚠ Files differ")
        result = compare_numeric_dfs(czid_sorted, seqtoid_sorted, id_cols=['sample_name'])
        print(f"    → {result}")

    # Check missing/extra samples
    for name, df in [('czid', czid_df), ('seqtoid', seqtoid_df)]:
        actual = set(df['sample_name'].astype(str).str.strip())
        missing = set(EXPECTED_SAMPLES) - actual
        extra   = actual - set(EXPECTED_SAMPLES)
        if missing:
            print(f"  ⚠ Missing in {name}: {sorted(missing)}")
        if extra:
            print(f"  ⚠ Extra in {name}:   {sorted(extra)}")


# ────────────────────────────────────────────────────────────────
# Step 2: Combined AMR Results
# ────────────────────────────────────────────────────────────────

def compare_combined_amr_results():
    file_name = 'combined_amr_results.csv'
    czid_path = os.path.join(CZID_DIR, file_name)
    seqtoid_path = os.path.join(SEQTOID_DIR, file_name)

    print("\n=== Step 2: Combined AMR Results Comparison ===")

    if not all(os.path.exists(p) for p in [czid_path, seqtoid_path]):
        print(f"✗ One or both {file_name} missing.")
        return

    czid_df    = pd.read_csv(czid_path)
    seqtoid_df = pd.read_csv(seqtoid_path)

    # Sort by the most stable composite key
    sort_cols = ['sample_name', 'gene_name']
    czid_sorted    = czid_df.sort_values(sort_cols).reset_index(drop=True)
    seqtoid_sorted = seqtoid_df.sort_values(sort_cols).reset_index(drop=True)

    print(f"Rows → czid: {len(czid_sorted):,} | seqtoid: {len(seqtoid_sorted):,}")

    # Quick set-of-genes-per-sample check
    czid_pairs = set(zip(czid_sorted['sample_name'], czid_sorted['gene_name']))
    seqtoid_pairs = set(zip(seqtoid_sorted['sample_name'], seqtoid_sorted['gene_name']))

    if czid_pairs != seqtoid_pairs:
        print("  ⚠ Detected gene-sample pairs differ!")
        if extra_c := czid_pairs - seqtoid_pairs:
            print(f"    Extra in czid: {len(extra_c)} pairs")
        if extra_s := seqtoid_pairs - czid_pairs:
            print(f"    Extra in seqtoid: {len(extra_s)} pairs")

    if czid_sorted.equals(seqtoid_sorted):
        print("  ✓ Files are identical (strict equality after sorting)")
    else:
        print("  ⚠ Strict comparison failed → numeric tolerance check")
        result = compare_numeric_dfs(
            czid_sorted, seqtoid_sorted,
            id_cols=sort_cols
        )
        print(f"    → {result}")

        # Highlight key AMR metrics if they show non-trivial differences
        key_metrics = ['rpm', 'dpm', 'read_coverage_depth', 'contig_coverage_breadth',
                       'contig_percent_id', 'num_reads', 'num_contigs']
        available = [c for c in key_metrics if c in czid_sorted.columns]
        for col in available:
            vals1 = czid_sorted[col].fillna(0).values
            vals2 = seqtoid_sorted[col].fillna(0).values
            cat = numeric_diff(vals1, vals2)
            if "warning" in cat or "significant" in cat:
                print(f"      {col:22} → {cat}")


# ────────────────────────────────────────────────────────────────
# Main
# ────────────────────────────────────────────────────────────────

def main():
    print("CZID AMR Pipeline Comparison (czid vs seqtoid)")
    print(f"  Focused on: sample_metadata.csv + combined_amr_results.csv")
    print(f"  Expected samples: {', '.join(EXPECTED_SAMPLES)}\n")

    compare_metadata()
    compare_combined_amr_results()

    print("\nComparison finished.")


if __name__ == '__main__':
    main()