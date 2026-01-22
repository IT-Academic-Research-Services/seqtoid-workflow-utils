#!/usr/bin/env python3
"""
CZID Consensus Genome Pipeline Comparison (czid vs seqtoid)

Compares outputs step-by-step between two pipeline runs.
Currently implemented:
  Step 1: sample_metadata.csv
  Step 2: consensus_genome_overview.csv
  Step 3: consensus FASTA files (~*_consensus.fa*)

Run this script from the directory that contains:
    czid/
    seqtoid/
"""

import os
import glob
import hashlib
import subprocess
import sys
import pandas as pd
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

METADATA_FILE           = 'sample_metadata.csv'
OVERVIEW_FILE           = 'consensus_genome_overview.csv'
CONSENSUS_PATTERN       = "*_consensus.fa*"   # catches .fa and .fasta


# ────────────────────────────────────────────────────────────────
# Shared helper functions
# ────────────────────────────────────────────────────────────────

def file_sha256(filepath: str) -> str:
    """Compute SHA-256 hash of a file in chunks."""
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

    sort_key = 'sample_name'
    if sort_key in czid_df.columns and sort_key in seqtoid_df.columns:
        czid_sorted    = czid_df.sort_values(sort_key).reset_index(drop=True)
        seqtoid_sorted = seqtoid_df.sort_values(sort_key).reset_index(drop=True)
    else:
        czid_sorted    = czid_df.reset_index(drop=True)
        seqtoid_sorted = seqtoid_df.reset_index(drop=True)

    if czid_sorted.equals(seqtoid_sorted):
        print("  ✓ Files are identical (strict equality after sorting)")
    else:
        print("  ⚠ Files differ")
        result = compare_numeric_dfs(czid_sorted, seqtoid_sorted, id_cols=[sort_key] if sort_key in czid_sorted else [])
        print(f"    → {result}")

        if set(czid_sorted.columns) != set(seqtoid_sorted.columns):
            print("      Column sets differ")
        if len(czid_sorted) != len(seqtoid_sorted):
            print(f"      Row count: czid={len(czid_sorted)}, seqtoid={len(seqtoid_sorted)}")

    for name, df in [('czid', czid_df), ('seqtoid', seqtoid_df)]:
        actual = set(df.get('sample_name', pd.Series()).astype(str).str.strip())
        expected = set(EXPECTED_SAMPLES)
        missing = expected - actual
        extra   = actual   - expected
        if missing:
            print(f"  ⚠ Missing in {name}: {', '.join(sorted(missing))}")
        if extra:
            print(f"  ⚠ Extra in {name}:   {', '.join(sorted(extra))}")


# ────────────────────────────────────────────────────────────────
# Step 2: Consensus Genome Overview
# ────────────────────────────────────────────────────────────────

def compare_consensus_overview():
    print("\n=== Step 2: Consensus Genome Overview Comparison ===")

    czid_path    = os.path.join(CZID_DIR,    OVERVIEW_FILE)
    seqtoid_path = os.path.join(SEQTOID_DIR, OVERVIEW_FILE)

    if not all(os.path.exists(p) for p in [czid_path, seqtoid_path]):
        print("✗ One or both files missing.")
        return

    try:
        czid_df    = pd.read_csv(czid_path)
        seqtoid_df = pd.read_csv(seqtoid_path)
    except Exception as e:
        print(f"✗ Failed to read one or both files: {e}")
        return

    sort_key = 'Sample Name'
    if sort_key in czid_df.columns and sort_key in seqtoid_df.columns:
        czid_sorted    = czid_df.sort_values(sort_key).reset_index(drop=True)
        seqtoid_sorted = seqtoid_df.sort_values(sort_key).reset_index(drop=True)
    else:
        czid_sorted    = czid_df.reset_index(drop=True)
        seqtoid_sorted = seqtoid_df.reset_index(drop=True)

    if czid_sorted.equals(seqtoid_sorted):
        print("  ✓ Files are identical (strict equality after sorting)")
    else:
        print("  ⚠ Files differ")
        result = compare_numeric_dfs(
            czid_sorted, seqtoid_sorted,
            id_cols=[sort_key] if sort_key in czid_sorted else [],
            atol=0.005
        )
        print(f"    → {result}")

        if len(czid_sorted) != len(seqtoid_sorted):
            print(f"      Row count mismatch: czid={len(czid_sorted)}, seqtoid={len(seqtoid_sorted)}")

    # Sample coverage check (just like metadata)
    for name, df in [('czid', czid_df), ('seqtoid', seqtoid_df)]:
        actual = set(df.get('Sample Name', pd.Series()).astype(str).str.strip())
        expected = set(EXPECTED_SAMPLES)
        missing = expected - actual
        extra   = actual   - expected
        if missing:
            print(f"  ⚠ Missing samples in {name}: {', '.join(sorted(missing))}")
        if extra:
            print(f"  ⚠ Extra samples in {name}:   {', '.join(sorted(extra))}")


# ────────────────────────────────────────────────────────────────
# Step 3: Consensus FASTA files
# ────────────────────────────────────────────────────────────────

def compare_consensus_fasta():
    print("\n=== Step 3: Consensus FASTA Comparison ===")

    missing_czid    = []
    missing_seqtoid = []
    problematic     = []

    for sample in EXPECTED_SAMPLES:
        print(f"  → {sample} ... ", end="")

        czid_files    = glob.glob(os.path.join(CZID_DIR,    f"{sample}{CONSENSUS_PATTERN}"))
        seqtoid_files = glob.glob(os.path.join(SEQTOID_DIR, f"{sample}{CONSENSUS_PATTERN}"))

        if len(czid_files) != 1:
            print("✗", end=" ")
            if len(czid_files) == 0:
                print("missing in czid")
                missing_czid.append(sample)
            else:
                print(f"multiple in czid ({len(czid_files)})")
            problematic.append(sample)
            continue

        if len(seqtoid_files) != 1:
            print("✗", end=" ")
            if len(seqtoid_files) == 0:
                print("missing in seqtoid")
                missing_seqtoid.append(sample)
            else:
                print(f"multiple in seqtoid ({len(seqtoid_files)})")
            problematic.append(sample)
            continue

        czid_fa    = czid_files[0]
        seqtoid_fa = seqtoid_files[0]

        czid_hash    = file_sha256(czid_fa)
        seqtoid_hash = file_sha256(seqtoid_fa)

        if czid_hash == seqtoid_hash:
            print("✓ byte-for-byte identical")
            continue

        print("⚠ byte-for-byte differ → checking content ... ", end="")

        try:
            czid_lines    = sum(1 for _ in open(czid_fa))
            seqtoid_lines = sum(1 for _ in open(seqtoid_fa))
            if czid_lines != seqtoid_lines:
                print(f"different line count (czid={czid_lines}, seqtoid={seqtoid_lines})")
                problematic.append(sample)
                continue
        except Exception as e:
            print(f"line count error: {e}")
            problematic.append(sample)
            continue

        try:
            cmd = "seqkit seq -s -i {} | sort | sha256sum"
            czid_seq   = subprocess.check_output(cmd.format(czid_fa),    shell=True, text=True).split()[0]
            seqtoid_seq = subprocess.check_output(cmd.format(seqtoid_fa), shell=True, text=True).split()[0]

            if czid_seq == seqtoid_seq:
                print("sequences match (order & formatting ignored)")
            else:
                print("✗ sequences DIFFER")
                problematic.append(sample)

        except FileNotFoundError:
            print("seqkit not found → skipping deep check")
            problematic.append(sample)
        except subprocess.CalledProcessError as e:
            print(f"seqkit failed: {e}")
            problematic.append(sample)

    if missing_czid or missing_seqtoid:
        if missing_czid:
            print(f"  Missing FASTA in czid:    {', '.join(sorted(missing_czid))}")
        if missing_seqtoid:
            print(f"  Missing FASTA in seqtoid: {', '.join(sorted(missing_seqtoid))}")

    if problematic:
        print(f"  Samples with FASTA issues: {', '.join(sorted(set(problematic)))}")


# ────────────────────────────────────────────────────────────────
# Main execution
# ────────────────────────────────────────────────────────────────

def main():
    print("CZID Consensus Genome Pipeline Comparison")
    print(f"  {CZID_DIR}/  vs  {SEQTOID_DIR}/")
    print(f"  Samples: {', '.join(EXPECTED_SAMPLES)}\n")

    compare_metadata()
    compare_consensus_overview()
    compare_consensus_fasta()

    print("\nComparison finished (steps 1–3).")


if __name__ == '__main__':
    main()