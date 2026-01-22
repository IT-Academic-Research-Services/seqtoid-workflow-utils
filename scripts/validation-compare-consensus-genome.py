#!/usr/bin/env python3
"""
CZID Consensus Genome Pipeline Comparison (czid vs seqtoid)

Currently implemented:
  Step 1: sample_metadata.csv
  Step 2: consensus_genome_overview.csv
  Step 3: per-sample assembly/QC report (report.tsv inside sample subfolder)
  Step 4: consensus FASTA files (~*_consensus.fa*)
  Step 5: per-sample samtools depth (samtools_depth.txt inside sample subfolder)

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
from typing import List, Optional

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
REPORT_FILE       = 'report.tsv'
DEPTH_FILE        = 'samtools_depth.txt'
CONSENSUS_PATTERN = "*_consensus.fa*"


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


def load_depth(path: str) -> Optional[pd.Series]:
    """Load samtools_depth.txt as a numeric Series (one depth per line)"""
    try:
        return pd.read_csv(path, header=None, dtype=float)[0]
    except Exception as e:
        print(f"  Failed to read depth file {path}: {e}", file=sys.stderr)
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
    czid_sorted    = czid_df.sort_values(sort_key).reset_index(drop=True)    if sort_key in czid_df.columns else czid_df.reset_index(drop=True)
    seqtoid_sorted = seqtoid_df.sort_values(sort_key).reset_index(drop=True) if sort_key in seqtoid_df.columns else seqtoid_df.reset_index(drop=True)

    if czid_sorted.equals(seqtoid_sorted):
        print("  ✓ Files are identical (strict equality after sorting)")
    else:
        print("  ⚠ Files differ")
        result = compare_numeric_dfs(czid_sorted, seqtoid_sorted, atol=0.005)
        print(f"    → {result}")


# ────────────────────────────────────────────────────────────────
# Step 3: Per-sample report.tsv
# ────────────────────────────────────────────────────────────────

def compare_quast_reports():
    print("\n=== Step 3: Per-sample assembly/QC report (report.tsv) Comparison ===")
    problematic = []

    for sample in EXPECTED_SAMPLES:
        print(f"  → {sample} ... ", end="")

        czid_subdirs    = glob.glob(os.path.join(CZID_DIR,    f"{sample}_*"))
        seqtoid_subdirs = glob.glob(os.path.join(SEQTOID_DIR, f"{sample}_*"))

        if len(czid_subdirs) != 1 or len(seqtoid_subdirs) != 1:
            print("✗ subfolder issue")
            problematic.append(sample)
            continue

        czid_report    = os.path.join(czid_subdirs[0], REPORT_FILE)
        seqtoid_report = os.path.join(seqtoid_subdirs[0], REPORT_FILE)

        if not os.path.isfile(czid_report) or not os.path.isfile(seqtoid_report):
            print("✗ report.tsv missing")
            problematic.append(sample)
            continue

        czid_df    = pd.read_csv(czid_report, sep='\t', header=None, names=['Metric', 'Value'], dtype=str)
        seqtoid_df = pd.read_csv(seqtoid_report, sep='\t', header=None, names=['Metric', 'Value'], dtype=str)

        czid_df['Value_num'] = pd.to_numeric(czid_df['Value'], errors='coerce')
        seqtoid_df['Value_num'] = pd.to_numeric(seqtoid_df['Value'], errors='coerce')

        merged = czid_df.merge(seqtoid_df, on='Metric', suffixes=('_czid', '_seqtoid'), how='outer')
        num_df = pd.DataFrame({
            'czid': merged['Value_num_czid'].dropna(),
            'seqtoid': merged['Value_num_seqtoid'].dropna()
        }).reset_index(drop=True)

        if len(num_df) == 0:
            print("no numeric values")
            continue

        result = compare_numeric_dfs(num_df, num_df, id_cols=[], atol=0.005)  # dummy second df to reuse function
        # Better: compare the two value columns
        comp_df = pd.DataFrame({'czid': merged['Value_num_czid'].fillna(0), 'seqtoid': merged['Value_num_seqtoid'].fillna(0)})
        result = compare_numeric_dfs(comp_df, comp_df, id_cols=[], atol=0.005)  # reuse logic

        if "equivalent" in result or "identical" in result:
            print("✓ equivalent within tolerance")
        elif "warning" in result:
            print(f"⚠ minor differences")
            problematic.append(sample)
        else:
            print(f"✗ {result}")
            problematic.append(sample)

    if problematic:
        print(f"\n  Samples with report issues: {', '.join(sorted(set(problematic)))}")


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
            print("✗ missing/multiple")
            if len(czid_files) == 0: missing_czid.append(sample)
            if len(seqtoid_files) == 0: missing_seqtoid.append(sample)
            problematic.append(sample)
            continue

        czid_fa = czid_files[0]
        seqtoid_fa = seqtoid_files[0]

        if file_sha256(czid_fa) == file_sha256(seqtoid_fa):
            print("✓ byte-for-byte identical")
            continue

        print("⚠ byte-for-byte differ → content check ... ", end="")

        try:
            czid_lines = sum(1 for _ in open(czid_fa))
            seqtoid_lines = sum(1 for _ in open(seqtoid_fa))
            if czid_lines != seqtoid_lines:
                print("different line count")
                problematic.append(sample)
                continue
        except Exception:
            print("line count error")
            problematic.append(sample)
            continue

        try:
            cmd = "seqkit seq -s -i {} | sort | sha256sum"
            h1 = subprocess.check_output(cmd.format(czid_fa), shell=True, text=True).split()[0]
            h2 = subprocess.check_output(cmd.format(seqtoid_fa), shell=True, text=True).split()[0]
            if h1 == h2:
                print("sequences match (order ignored)")
            else:
                print("✗ sequences DIFFER")
                problematic.append(sample)
        except FileNotFoundError:
            print("seqkit not found")
            problematic.append(sample)
        except Exception as e:
            print(f"seqkit error: {e}")
            problematic.append(sample)

    if missing_czid or missing_seqtoid:
        print(f"  Missing FASTA → czid: {missing_czid} | seqtoid: {missing_seqtoid}")


# ────────────────────────────────────────────────────────────────
# Step 5: samtools_depth.txt (per-position depth)
# ────────────────────────────────────────────────────────────────

def compare_samtools_depth():
    print("\n=== Step 5: samtools depth (samtools_depth.txt) Comparison ===")
    problematic = []

    for sample in EXPECTED_SAMPLES:
        print(f"  → {sample} ... ", end="")

        czid_subdirs    = glob.glob(os.path.join(CZID_DIR,    f"{sample}_*"))
        seqtoid_subdirs = glob.glob(os.path.join(SEQTOID_DIR, f"{sample}_*"))

        if len(czid_subdirs) != 1 or len(seqtoid_subdirs) != 1:
            print("✗ subfolder issue")
            problematic.append(sample)
            continue

        czid_depth_path    = os.path.join(czid_subdirs[0],    DEPTH_FILE)
        seqtoid_depth_path = os.path.join(seqtoid_subdirs[0], DEPTH_FILE)

        if not os.path.isfile(czid_depth_path) or not os.path.isfile(seqtoid_depth_path):
            print("✗ depth file missing")
            problematic.append(sample)
            continue

        czid_depth    = load_depth(czid_depth_path)
        seqtoid_depth = load_depth(seqtoid_depth_path)

        if czid_depth is None or seqtoid_depth is None:
            print("failed to load")
            problematic.append(sample)
            continue

        if len(czid_depth) != len(seqtoid_depth):
            print(f"different length (czid={len(czid_depth)}, seqtoid={len(seqtoid_depth)})")
            problematic.append(sample)
            continue

        comp_df = pd.DataFrame({
            'czid': czid_depth.values,
            'seqtoid': seqtoid_depth.values
        })

        result = compare_numeric_dfs(comp_df, comp_df, id_cols=[], atol=0.005)

        if "equivalent" in result or "identical" in result:
            print("✓ depths equivalent within tolerance")
        elif "warning" in result:
            print(f"⚠ minor differences → {result}")
            problematic.append(sample)
        else:
            print(f"✗ {result}")
            problematic.append(sample)

    if problematic:
        print(f"\n  Samples with depth differences/issues: {', '.join(sorted(set(problematic)))}")


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
    compare_samtools_depth()

    print("\nComparison finished.")


if __name__ == '__main__':
    main()