#!/usr/bin/env python3
"""
CZID AMR Pipeline Comparison
Compares:
  1. sample_metadata.csv
  2. combined_amr_results.csv
  3. primary_AMR_report.tsv (per sample, inside final_reports/)
"""

import os
import glob
import pandas as pd
import numpy as np
import hashlib

# ────────────────────────────────────────────────────────────────
# Configuration
# ────────────────────────────────────────────────────────────────

CZID_DIR    = 'czid'
SEQTOID_DIR = 'seqtoid'

EXPECTED_SAMPLES = [
    'ERR11417004',
    'SRR10903401',
    'SRR12876565',
    'SRR13227003',
    'SRR13227004',
    'SRR13227005',
    'SRR15049352',
]

# ────────────────────────────────────────────────────────────────
# Helpers
# ────────────────────────────────────────────────────────────────

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
        id_cols=None,
        atol: float = 0.005
) -> str:
    if id_cols is None:
        id_cols = []
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
    file_path = 'sample_metadata.csv'
    czid_p = os.path.join(CZID_DIR, file_path)
    seq_p  = os.path.join(SEQTOID_DIR, file_path)

    if not all(os.path.exists(p) for p in [czid_p, seq_p]):
        print("✗ One or both sample_metadata.csv missing.")
        return

    czid_df = pd.read_csv(czid_p)
    seq_df  = pd.read_csv(seq_p)

    czid_s = czid_df.sort_values('sample_name').reset_index(drop=True)
    seq_s  = seq_df.sort_values('sample_name').reset_index(drop=True)

    if czid_s.equals(seq_s):
        print("  ✓ identical (after sorting)")
    else:
        print("  ⚠ differ → numeric check")
        print("    → " + compare_numeric_dfs(czid_s, seq_s, id_cols=['sample_name']))

    for name, df in [('czid', czid_df), ('seqtoid', seq_df)]:
        actual = set(df['sample_name'].astype(str).str.strip())
        miss = sorted(set(EXPECTED_SAMPLES) - actual)
        extra = sorted(actual - set(EXPECTED_SAMPLES))
        if miss:  print(f"  ⚠ Missing in {name}: {miss}")
        if extra: print(f"  ⚠ Extra in {name}:   {extra}")


# ────────────────────────────────────────────────────────────────
# Step 2: Combined AMR Results
# ────────────────────────────────────────────────────────────────

def compare_combined_amr_results():
    print("\n=== Step 2: combined_amr_results.csv Comparison ===")
    file_name = 'combined_amr_results.csv'
    czid_p = os.path.join(CZID_DIR, file_name)
    seq_p  = os.path.join(SEQTOID_DIR, file_name)

    if not all(os.path.exists(p) for p in [czid_p, seq_p]):
        print(f"✗ One or both {file_name} missing.")
        return

    czid_df = pd.read_csv(czid_p)
    seq_df  = pd.read_csv(seq_p)

    sort_cols = ['sample_name', 'gene_name']
    czid_s = czid_df.sort_values(sort_cols).reset_index(drop=True)
    seq_s  = seq_df.sort_values(sort_cols).reset_index(drop=True)

    print(f"Rows → czid: {len(czid_s):,} | seqtoid: {len(seq_s):,}")

    czid_pairs = set(zip(czid_s['sample_name'], czid_s['gene_name']))
    seq_pairs  = set(zip(seq_s['sample_name'], seq_s['gene_name']))

    if czid_pairs != seq_pairs:
        print("  ⚠ Gene-sample pairs differ!")
        if czid_pairs - seq_pairs:
            print(f"    Extra in czid: {len(czid_pairs - seq_pairs)}")
        if seq_pairs - czid_pairs:
            print(f"    Extra in seqtoid: {len(seq_pairs - czid_pairs)}")

    if czid_s.equals(seq_s):
        print("  ✓ identical (after sorting)")
    else:
        print("  ⚠ differ → numeric check")
        result = compare_numeric_dfs(czid_s, seq_s, id_cols=sort_cols)
        print(f"    → {result}")

        # Highlight important columns
        key_cols = ['rpm', 'dpm', 'read_coverage_depth', 'contig_coverage_breadth',
                    'contig_percent_id', 'num_reads', 'num_contigs']
        for col in [c for c in key_cols if c in czid_s.columns]:
            cat = numeric_diff(
                czid_s[col].fillna(0).values,
                seq_s[col].fillna(0).values
            )
            if "warning" in cat or "significant" in cat:
                print(f"      {col:22} → {cat}")


# ────────────────────────────────────────────────────────────────
# Step 3: Per-sample primary_AMR_report.tsv
# ────────────────────────────────────────────────────────────────

def compare_primary_amr_reports():
    print("\n=== Step 3: Per-sample primary_AMR_report.tsv ===")
    missing_c = []
    missing_s = []

    for sample in EXPECTED_SAMPLES:
        print(f"  {sample} ... ", end="")

        cz_folders = glob.glob(os.path.join(CZID_DIR, f"{sample}_*"))
        seq_folders = glob.glob(os.path.join(SEQTOID_DIR, f"{sample}_*"))

        if len(cz_folders) != 1 or len(seq_folders) != 1:
            print("folder issue")
            if len(cz_folders) != 1: missing_c.append(sample)
            if len(seq_folders) != 1: missing_s.append(sample)
            continue

        cz_report = os.path.join(cz_folders[0], "final_reports", "primary_AMR_report.tsv")
        seq_report = os.path.join(seq_folders[0], "final_reports", "primary_AMR_report.tsv")

        if not os.path.isfile(cz_report) or not os.path.isfile(seq_report):
            print("missing report")
            if not os.path.isfile(cz_report): missing_c.append(sample)
            if not os.path.isfile(seq_report): missing_s.append(sample)
            continue

        try:
            cz_df = pd.read_csv(cz_report, sep='\t')
            seq_df = pd.read_csv(seq_report, sep='\t')

            sort_col = 'gene_name'
            cz_s = cz_df.sort_values(sort_col).reset_index(drop=True)
            seq_s = seq_df.sort_values(sort_col).reset_index(drop=True)

            if cz_s.equals(seq_s):
                print("identical")
            else:
                print("differ → numeric check")
                result = compare_numeric_dfs(cz_s, seq_s, id_cols=[sort_col])
                print(f"      → {result}")
                if len(cz_s) != len(seq_s):
                    print(f"      rows: czid={len(cz_s)}, seqtoid={len(seq_s)}")

        except Exception as e:
            print(f"error: {e}")

    if missing_c:
        print(f"\n  Missing report in czid: {', '.join(missing_c)}")
    if missing_s:
        print(f"  Missing report in seqtoid: {', '.join(missing_s)}")


# ────────────────────────────────────────────────────────────────
# Main
# ────────────────────────────────────────────────────────────────

def main():
    print("CZID AMR Pipeline Comparison (czid vs seqtoid)")
    print("  Files checked:")
    print("    • sample_metadata.csv")
    print("    • combined_amr_results.csv")
    print("    • primary_AMR_report.tsv (per sample)")
    print(f"  Samples: {', '.join(EXPECTED_SAMPLES)}\n")

    compare_metadata()
    compare_combined_amr_results()
    compare_primary_amr_reports()

    print("\nComparison finished.")


if __name__ == '__main__':
    main()