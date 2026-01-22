#!/usr/bin/env python3
"""
CZID Consensus Genome Pipeline Comparison (czid vs seqtoid)

Currently implemented:
  Step 1: sample_metadata.csv
  Step 2: consensus_genome_overview.csv
  Step 3: per-sample assembly/QC report (report.tsv inside SRR... subfolder)
  Step 4: consensus FASTA files (~*_consensus.fa*)

Run from directory containing:
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
from typing import List, Optional, Tuple

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

METADATA_FILE     = 'sample_metadata.csv'
OVERVIEW_FILE     = 'consensus_genome_overview.csv'
REPORT_FILE       = 'report.tsv'               # inside sample subfolder
CONSENSUS_PATTERN = "*_consensus.fa*"          # top-level or inside subfolder


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


def load_key_value_tsv(path: str) -> Optional[pd.DataFrame]:
    """Load QUAST-style report.tsv as key-value DataFrame"""
    try:
        df = pd.read_csv(path, sep='\t', header=None, names=['Metric', 'Value'], dtype=str)
        df['Value_num'] = pd.to_numeric(df['Value'], errors='coerce')
        return df
    except Exception as e:
        print(f"  Failed to read {path}: {e}", file=sys.stderr)
        return None


# ────────────────────────────────────────────────────────────────
# Step 1: Sample Metadata
# ────────────────────────────────────────────────────────────────

def compare_metadata():
    print("\n=== Step 1: Sample Metadata Comparison ===")
    czid_path    = os.path.join(CZID_DIR,    METADATA_FILE)
    seqtoid_path = os.path.join(SEQTOID_DIR, METADATA_FILE)

    if not all(os.path.exists(p) for p in [czid_path, seqtoid_path]):
        print("✗ One or both files missing.")
        return

    czid_df    = pd.read_csv(czid_path, dtype=str)
    seqtoid_df = pd.read_csv(seqtoid_path, dtype=str)

    sort_key = 'sample_name'
    czid_sorted    = czid_df.sort_values(sort_key).reset_index(drop=True)    if sort_key in czid_df.columns else czid_df.reset_index(drop=True)
    seqtoid_sorted = seqtoid_df.sort_values(sort_key).reset_index(drop=True) if sort_key in seqtoid_df.columns else seqtoid_df.reset_index(drop=True)

    if czid_sorted.equals(seqtoid_sorted):
        print("  ✓ Files are identical (strict equality after sorting)")
    else:
        print("  ⚠ Files differ")
        result = compare_numeric_dfs(czid_sorted, seqtoid_sorted, atol=0.005)
        print(f"    → {result}")

    for name, df in [('czid', czid_df), ('seqtoid', seqtoid_df)]:
        actual = set(df.get('sample_name', pd.Series()).astype(str).str.strip())
        missing = set(EXPECTED_SAMPLES) - actual
        extra   = actual - set(EXPECTED_SAMPLES)
        if missing: print(f"  ⚠ Missing in {name}: {sorted(missing)}")
        if extra:   print(f"  ⚠ Extra in {name}:   {sorted(extra)}")


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

    czid_df    = pd.read_csv(czid_path)
    seqtoid_df = pd.read_csv(seqtoid_path)

    sort_key = 'Sample Name'
    czid_sorted    = czid_df.sort_values(sort_key).reset_index(drop=True)    if sort_key in czid_df else czid_df.reset_index(drop=True)
    seqtoid_sorted = seqtoid_df.sort_values(sort_key).reset_index(drop=True) if sort_key in seqtoid_df else seqtoid_df.reset_index(drop=True)

    if czid_sorted.equals(seqtoid_sorted):
        print("  ✓ Files are identical (strict equality after sorting)")
    else:
        print("  ⚠ Files differ")
        result = compare_numeric_dfs(czid_sorted, seqtoid_sorted, atol=0.005)
        print(f"    → {result}")


# ────────────────────────────────────────────────────────────────
# Step 3: Per-sample QUAST-style report.tsv
# ────────────────────────────────────────────────────────────────

def compare_quast_reports():
    print("\n=== Step 3: Per-sample assembly/QC report (report.tsv) Comparison ===")

    problematic = []

    for sample in EXPECTED_SAMPLES:
        print(f"  → {sample} ... ", end="")

        # Look for subfolder like SRR34692681_*
        czid_subdirs    = glob.glob(os.path.join(CZID_DIR,    f"{sample}_*"))
        seqtoid_subdirs = glob.glob(os.path.join(SEQTOID_DIR, f"{sample}_*"))

        if len(czid_subdirs) != 1 or len(seqtoid_subdirs) != 1:
            print("✗ subfolder issue")
            problematic.append(sample)
            continue

        czid_report    = os.path.join(czid_subdirs[0],    REPORT_FILE)
        seqtoid_report = os.path.join(seqtoid_subdirs[0], REPORT_FILE)

        if not os.path.isfile(czid_report) or not os.path.isfile(seqtoid_report):
            print("✗ report.tsv missing in one or both")
            problematic.append(sample)
            continue

        czid_df    = load_key_value_tsv(czid_report)
        seqtoid_df = load_key_value_tsv(seqtoid_report)

        if czid_df is None or seqtoid_df is None:
            print("failed to parse")
            problematic.append(sample)
            continue

        # Merge on Metric to allow side-by-side comparison
        merged = czid_df.merge(seqtoid_df, on='Metric', suffixes=('_czid', '_seqtoid'), how='outer')

        # Numeric comparison
        numeric_cols = ['Value_num_czid', 'Value_num_seqtoid']
        num_merged = merged[numeric_cols].dropna(how='any')

        if len(num_merged) == 0:
            print("no comparable numeric values")
            continue

        result = compare_numeric_dfs(
            num_merged.rename(columns={'Value_num_czid': 'czid', 'Value_num_seqtoid': 'seqtoid'}),
            id_cols=[],
            atol=0.005
        )

        if "equivalent" in result or "identical" in result:
            print("✓ equivalent within tolerance")
        elif "warning" in result:
            print(f"⚠ minor differences → {result}")
            problematic.append(sample)
        else:
            print(f"✗ {result}")
            problematic.append(sample)

    if problematic:
        print(f"\n  Samples with report differences / issues: {', '.join(sorted(set(problematic)))}")


# ────────────────────────────────────────────────────────────────
# Step 4: Consensus FASTA
# ────────────────────────────────────────────────────────────────

def compare_consensus_fasta():
    print("\n=== Step 4: Consensus FASTA Comparison ===")
    missing_czid = []
    missing_seqtoid = []
    problematic = []

    for sample in EXPECTED_SAMPLES:
        print(f"  → {sample} ... ", end="")

        czid_files = glob.glob(os.path.join(CZID_DIR, f"{sample}{CONSENSUS_PATTERN}"))
        seqtoid_files = glob.glob(os.path.join(SEQTOID_DIR, f"{sample}{CONSENSUS_PATTERN}"))

        if len(czid_files) != 1 or len(seqtoid_files) != 1:
            print("✗ file count mismatch / missing")
            if len(czid_files) == 0: missing_czid.append(sample)
            if len(seqtoid_files) == 0: missing_seqtoid.append(sample)
            problematic.append(sample)
            continue

        czid_fa = czid_files[0]
        seqtoid_fa = seqtoid_files[0]

        if file_sha256(czid_fa) == file_sha256(seqtoid_fa):
            print("✓ byte-for-byte identical")
            continue

        print("⚠ differ byte-for-byte → content check ... ", end="")

        try:
            czid_lines = sum(1 for _ in open(czid_fa))
            seqtoid_lines = sum(1 for _ in open(seqtoid_fa))
            if czid_lines != seqtoid_lines:
                print("different line count")
                problematic.append(sample)
                continue
        except:
            print("line count error")
            problematic.append(sample)
            continue

        try:
            cmd = "seqkit seq -s -i {} | sort | sha256sum"
            h1 = subprocess.check_output(cmd.format(czid_fa),    shell=True, text=True).split()[0]
            h2 = subprocess.check_output(cmd.format(seqtoid_fa), shell=True, text=True).split()[0]
            print("sequences match (order ignored)" if h1 == h2 else "✗ sequences DIFFER")
        except FileNotFoundError:
            print("seqkit not found")
            problematic.append(sample)
        except Exception as e:
            print(f"seqkit error: {e}")
            problematic.append(sample)

    if missing_czid or missing_seqtoid:
        print(f"  Missing FASTA → czid: {missing_czid} | seqtoid: {missing_seqtoid}")


# ────────────────────────────────────────────────────────────────
# Main
# ────────────────────────────────────────────────────────────────

def main():
    print("CZID Consensus Genome Pipeline Comparison")
    print(f"  {CZID_DIR}/ vs {SEQTOID_DIR}/")
    print(f"  Samples: {', '.join(EXPECTED_SAMPLES)}\n")

    compare_metadata()
    compare_consensus_overview()
    compare_quast_reports()
    compare_consensus_fasta()

    print("\nComparison finished.")


if __name__ == '__main__':
    main()