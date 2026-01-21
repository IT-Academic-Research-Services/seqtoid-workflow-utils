import pandas as pd
import os
import glob
import numpy as np
import hashlib
import subprocess  # for future FASTQ/FASTA steps

# Define paths to the directories
CZID_DIR = 'czid'
SEQTOID_DIR = 'seqtoid'

# List of expected sample IDs
EXPECTED_SAMPLES = [
    'SRR18291896',
    'SRR15049352',
    'SRR12048509'
]

# ────────────────────────────────────────────────────────────────
# Helper functions
# ────────────────────────────────────────────────────────────────

def file_sha256(filepath):
    """Compute SHA-256 hash of a file in chunks."""
    sha256 = hashlib.sha256()
    with open(filepath, 'rb') as f:
        while chunk := f.read(8192 * 1024):
            sha256.update(chunk)
    return sha256.hexdigest()

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
        id_cols: list,
        atol: float = 0.005
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
# Step 1: Sample Metadata Comparison
# ────────────────────────────────────────────────────────────────

def compare_metadata():
    metadata_file = 'sample_metadata.csv'
    czid_path = os.path.join(CZID_DIR, metadata_file)
    seqtoid_path = os.path.join(SEQTOID_DIR, metadata_file)

    print("=== Step 1: Sample Metadata Comparison ===")

    if not os.path.exists(czid_path):
        print(f"✗ Missing in czid: {czid_path}")
        return
    if not os.path.exists(seqtoid_path):
        print(f"✗ Missing in seqtoid: {seqtoid_path}")
        return

    czid_df = pd.read_csv(czid_path)
    seqtoid_df = pd.read_csv(seqtoid_path)

    czid_sorted = czid_df.sort_values('sample_name').reset_index(drop=True)
    seqtoid_sorted = seqtoid_df.sort_values('sample_name').reset_index(drop=True)

    if czid_sorted.equals(seqtoid_sorted):
        print("  ✓ Files are identical (strict equality after sorting)")
    else:
        print("  ⚠ Files differ")
        if list(czid_sorted.columns) != list(seqtoid_sorted.columns):
            print("    Columns differ!")
            print("      czid:   ", list(czid_sorted.columns))
            print("      seqtoid:", list(seqtoid_sorted.columns))
        try:
            result = compare_numeric_dfs(czid_sorted, seqtoid_sorted, id_cols=['sample_name'])
            print(f"    → {result}")
        except Exception as e:
            print(f"    → Tolerant numeric check error: {e}")

    for name, df in [('czid', czid_df), ('seqtoid', seqtoid_df)]:
        actual = set(df['sample_name'].astype(str).str.strip())
        expected = set(EXPECTED_SAMPLES)
        missing = expected - actual
        extra = actual - expected
        if missing:
            print(f"  ⚠ Missing in {name}: {sorted(missing)}")
        if extra:
            print(f"  ⚠ Extra in {name}: {sorted(extra)}")


# ────────────────────────────────────────────────────────────────
# Step 2: Per-sample Taxon Reports Comparison
# ────────────────────────────────────────────────────────────────

def compare_taxon_reports():
    print("\n=== Step 2: Sample Taxon Reports Comparison ===")
    missing_czid = []
    missing_seqtoid = []

    for sample in EXPECTED_SAMPLES:
        czid_files = glob.glob(os.path.join(CZID_DIR, f"{sample}_*_taxon_report.csv"))
        seqtoid_files = glob.glob(os.path.join(SEQTOID_DIR, f"{sample}_*_taxon_report.csv"))

        if len(czid_files) != 1:
            missing_czid.append(sample)
            print(f"  ✗ {sample}: missing or multiple files in czid")
            continue
        if len(seqtoid_files) != 1:
            missing_seqtoid.append(sample)
            print(f"  ✗ {sample}: missing or multiple files in seqtoid")
            continue

        czid_path = czid_files[0]
        seqtoid_path = seqtoid_files[0]

        print(f"  → {sample} ...", end=" ")

        try:
            czid_df = pd.read_csv(czid_path, dtype=str)
            seqtoid_df = pd.read_csv(seqtoid_path, dtype=str)

            # Convert tax_id for sorting
            czid_df['tax_id'] = pd.to_numeric(czid_df['tax_id'], errors='coerce')
            seqtoid_df['tax_id'] = pd.to_numeric(seqtoid_df['tax_id'], errors='coerce')

            czid_sorted = czid_df.sort_values('tax_id').reset_index(drop=True)
            seqtoid_sorted = seqtoid_df.sort_values('tax_id').reset_index(drop=True)

            if czid_sorted.equals(seqtoid_sorted):
                print("identical")
            else:
                print("strict fail → tolerant numeric check")
                result = compare_numeric_dfs(czid_sorted, seqtoid_sorted, id_cols=['tax_id'])
                print(f"      → {result}")

                # Optional: show row count difference
                if len(czid_sorted) != len(seqtoid_sorted):
                    print(f"      Row count: czid={len(czid_sorted)}, seqtoid={len(seqtoid_sorted)}")

        except Exception as e:
            print(f"failed to compare: {e}")

    if missing_czid:
        print(f"  Missing taxon_report files in czid for: {', '.join(missing_czid)}")
    if missing_seqtoid:
        print(f"  Missing taxon_report files in seqtoid for: {', '.join(missing_seqtoid)}")


# ────────────────────────────────────────────────────────────────
# Placeholders for future steps (expand as you identify files)
# ────────────────────────────────────────────────────────────────

def compare_step_3():
    print("\n=== Step 3: [TODO: Assembly / QC Statistics] ===")
    print("  Placeholder: Compare assembly stats, N50, total contigs, read counts, etc.")

def compare_step_4():
    print("\n=== Step 4: [TODO: Combined Taxon Results / RPM tables] ===")
    print("  Placeholder: Compare combined_sample_taxon_results_NT.rpm.csv or equivalent")

def compare_nonhost_fastqs():
    print("\n=== Step X: Non-host reads (FASTQ) ===")
    print("  Placeholder: Implement when ready (byte hash + seqkit sequence hash)")

def compare_nonhost_contigs():
    print("\n=== Step Y: Non-host contigs (FASTA) ===")
    print("  Placeholder: Implement when ready")


def main():
    print("Starting CZID Long Reads Pipeline Comparison (czid vs seqtoid)\n")
    compare_metadata()
    compare_taxon_reports()
    compare_step_3()
    compare_step_4()
    # compare_nonhost_fastqs()
    # compare_nonhost_contigs()
    print("\nComparison complete. Add next step functions as needed.")


if __name__ == '__main__':
    main()