import pandas as pd
import os
import sys
import numpy as np
from typing import List

# ────────────────────────────────────────────────────────────────
# Configuration
# ────────────────────────────────────────────────────────────────

CZID_DIR    = 'czid'
SEQTOID_DIR = 'seqtoid'

EXPECTED_SAMPLES: List[str] = [
    'SRR11454628',
    'SRR34692683',
    'SRR34692681'
]

METADATA_FILE = 'sample_metadata.csv'

# ────────────────────────────────────────────────────────────────
# Helpers  (copied/adapted from long-reads version for consistency)
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
        id_cols: list,
        atol: float = 0.005,
) -> str:
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
    sig_cols    = [col for col, cat in results.items() if "significant" in cat]

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

    czid_path    = os.path.join(CZID_DIR,    METADATA_FILE)
    seqtoid_path = os.path.join(SEQTOID_DIR, METADATA_FILE)

    if not all(os.path.exists(p) for p in [czid_path, seqtoid_path]):
        print("✗ One or both files missing.")
        for p in [czid_path, seqtoid_path]:
            if not os.path.exists(p):
                print(f"  Missing: {p}")
        return

    try:
        czid_df    = pd.read_csv(czid_path,    dtype=str)
        seqtoid_df = pd.read_csv(seqtoid_path, dtype=str)
    except Exception as e:
        print(f"✗ Failed to read one or both files: {e}")
        return

    # Sort by sample_name (most common key in CZID metadata)
    sort_key = 'sample_name'
    if sort_key not in czid_df.columns or sort_key not in seqtoid_df.columns:
        print(f"⚠ Cannot sort — missing '{sort_key}' column in one or both files")
        czid_sorted    = czid_df.reset_index(drop=True)
        seqtoid_sorted = seqtoid_df.reset_index(drop=True)
    else:
        czid_sorted    = czid_df.sort_values(sort_key).reset_index(drop=True)
        seqtoid_sorted = seqtoid_df.sort_values(sort_key).reset_index(drop=True)

    if czid_sorted.equals(seqtoid_sorted):
        print("  ✓ Files are identical (strict equality after sorting)")
    else:
        print("  ⚠ Files differ")
        # Try tolerant numeric comparison anyway (for consistency with other steps)
        result = compare_numeric_dfs(czid_sorted, seqtoid_sorted, id_cols=[sort_key])
        print(f"    → {result}")

        # Also report obvious structural differences
        if list(czid_sorted.columns) != list(seqtoid_sorted.columns):
            print("      Column mismatch:")
            print(f"        czid:    {list(czid_sorted.columns)}")
            print(f"        seqtoid: {list(seqtoid_sorted.columns)}")
        if len(czid_sorted) != len(seqtoid_sorted):
            print(f"      Row count mismatch: czid={len(czid_sorted)}, seqtoid={len(seqtoid_sorted)}")

    # Sample coverage check
    for name, df in [('czid', czid_df), ('seqtoid', seqtoid_df)]:
        actual = set(df.get('sample_name', pd.Series()).astype(str).str.strip())
        expected = set(EXPECTED_SAMPLES)

        missing = expected - actual
        extra   = actual   - expected

        if missing:
            print(f"  ⚠ Missing samples in {name}: {', '.join(sorted(missing))}")
        if extra:
            print(f"  ⚠ Extra samples in {name}:   {', '.join(sorted(extra))}")


def main():
    print("CZID Consensus Genome Pipeline Comparison (czid vs seqtoid)")
    print(f"Expected samples: {', '.join(EXPECTED_SAMPLES)}\n")

    compare_metadata()

    # Future steps go here
    # print("\n=== Step 2: Consensus FASTA / assembly ===")
    # compare_consensus_fasta()

    print("\nComparison complete (so far).")


if __name__ == '__main__':
    main()