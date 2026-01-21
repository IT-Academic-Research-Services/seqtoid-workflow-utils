import pandas as pd
import os
import glob
import numpy as np
import hashlib
import subprocess

# ────────────────────────────────────────────────────────────────
# Paths & constants
# ────────────────────────────────────────────────────────────────

CZID_DIR = 'czid'
SEQTOID_DIR = 'seqtoid'

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
        wide_matrix: bool = False
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
    sig_cols = [col for col, cat in results.items() if "significant" in cat]

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
    # (unchanged – same as before)
    metadata_file = 'sample_metadata.csv'
    czid_path = os.path.join(CZID_DIR, metadata_file)
    seqtoid_path = os.path.join(SEQTOID_DIR, metadata_file)

    print("=== Step 1: Sample Metadata Comparison ===")
    if not all(os.path.exists(p) for p in [czid_path, seqtoid_path]):
        print("✗ One or both files missing.")
        return

    czid_df = pd.read_csv(czid_path)
    seqtoid_df = pd.read_csv(seqtoid_path)

    czid_sorted = czid_df.sort_values('sample_name').reset_index(drop=True)
    seqtoid_sorted = seqtoid_df.sort_values('sample_name').reset_index(drop=True)

    if czid_sorted.equals(seqtoid_sorted):
        print("  ✓ Files are identical (strict equality after sorting)")
    else:
        print("  ⚠ Files differ")
        result = compare_numeric_dfs(czid_sorted, seqtoid_sorted, id_cols=['sample_name'])
        print(f"    → {result}")

    for name, df in [('czid', czid_df), ('seqtoid', seqtoid_df)]:
        actual = set(df['sample_name'].astype(str).str.strip())
        expected = set(EXPECTED_SAMPLES)
        if missing := expected - actual:
            print(f"  ⚠ Missing in {name}: {sorted(missing)}")
        if extra := actual - expected:
            print(f"  ⚠ Extra in {name}: {sorted(extra)}")


# Step 2: Sample Overviews
# (unchanged)

def compare_sample_overviews():
    file_name = 'sample_overviews.csv'
    czid_path = os.path.join(CZID_DIR, file_name)
    seqtoid_path = os.path.join(SEQTOID_DIR, file_name)

    print("\n=== Step 2: Sample Overviews Comparison ===")
    if not all(os.path.exists(p) for p in [czid_path, seqtoid_path]):
        print(f"✗ One or both {file_name} missing.")
        return

    czid_df = pd.read_csv(czid_path, dtype=str)
    seqtoid_df = pd.read_csv(seqtoid_path, dtype=str)

    czid_sorted = czid_df.sort_values('sample_name').reset_index(drop=True)
    seqtoid_sorted = seqtoid_df.sort_values('sample_name').reset_index(drop=True)

    if czid_sorted.equals(seqtoid_sorted):
        print("  ✓ Files are identical (strict equality after sorting)")
    else:
        print("  ⚠ Strict comparison failed → checking numeric tolerance...")
        result = compare_numeric_dfs(czid_sorted, seqtoid_sorted, id_cols=['sample_name'])
        print("    → " + result)


# Step 3: Per-sample Taxon Reports
# (unchanged)

def compare_taxon_reports():
    print("\n=== Step 3: Per-sample Taxon Reports Comparison ===")
    missing_czid = []
    missing_seqtoid = []

    for sample in EXPECTED_SAMPLES:
        czid_files = glob.glob(os.path.join(CZID_DIR, f"{sample}_*_taxon_report.csv"))
        seqtoid_files = glob.glob(os.path.join(SEQTOID_DIR, f"{sample}_*_taxon_report.csv"))

        if len(czid_files) != 1:
            missing_czid.append(sample)
            print(f"  ✗ {sample}: missing/multiple taxon_report in czid")
            continue
        if len(seqtoid_files) != 1:
            missing_seqtoid.append(sample)
            print(f"  ✗ {sample}: missing/multiple taxon_report in seqtoid")
            continue

        czid_path = czid_files[0]
        seqtoid_path = seqtoid_files[0]
        print(f"  → {sample} ...", end=" ")

        try:
            czid_df = pd.read_csv(czid_path, dtype=str)
            seqtoid_df = pd.read_csv(seqtoid_path, dtype=str)

            czid_df['tax_id'] = pd.to_numeric(czid_df['tax_id'], errors='coerce')
            seqtoid_df['tax_id'] = pd.to_numeric(seqtoid_df['tax_id'], errors='coerce')

            czid_sorted = czid_df.sort_values('tax_id').reset_index(drop=True)
            seqtoid_sorted = seqtoid_df.sort_values('tax_id').reset_index(drop=True)

            if czid_sorted.equals(seqtoid_sorted):
                print("identical")
            else:
                print("strict fail → tolerant check")
                result = compare_numeric_dfs(czid_sorted, seqtoid_sorted, id_cols=['tax_id'])
                print(f"      → {result}")
                if len(czid_sorted) != len(seqtoid_sorted):
                    print(f"      Row count: czid={len(czid_sorted)}, seqtoid={len(seqtoid_sorted)}")
        except Exception as e:
            print(f"failed: {e}")

    if missing_czid:
        print(f"  Missing in czid: {', '.join(missing_czid)}")
    if missing_seqtoid:
        print(f"  Missing in seqtoid: {', '.join(missing_seqtoid)}")


# Step 4: Combined Taxon Results
# (unchanged)

def compare_combined_taxon_results():
    file_name = 'combined_sample_taxon_results_NT.bpm.csv'
    czid_path = os.path.join(CZID_DIR, file_name)
    seqtoid_path = os.path.join(SEQTOID_DIR, file_name)

    print("\n=== Step 4: Combined Sample Taxon Results (NT.bpm.csv) Comparison ===")

    if not all(os.path.exists(p) for p in [czid_path, seqtoid_path]):
        print(f"✗ One or both {file_name} missing.")
        return

    czid_df = pd.read_csv(czid_path)
    seqtoid_df = pd.read_csv(seqtoid_path)

    sort_col = 'Taxon Name'
    czid_sorted = czid_df.sort_values(sort_col).reset_index(drop=True)
    seqtoid_sorted = seqtoid_df.sort_values(sort_col).reset_index(drop=True)

    print(f"Taxa count - czid: {len(czid_sorted)}, seqtoid: {len(seqtoid_sorted)}")

    czid_taxa = set(czid_sorted[sort_col])
    seqtoid_taxa = set(seqtoid_sorted[sort_col])
    if czid_taxa != seqtoid_taxa:
        print("  ⚠ Taxon sets differ!")
        if extra := czid_taxa - seqtoid_taxa:
            print(f"    Extra in czid: {len(extra)} taxa")
        if extra := seqtoid_taxa - czid_taxa:
            print(f"    Extra in seqtoid: {len(extra)} taxa")

    if czid_sorted.equals(seqtoid_sorted):
        print("  ✓ Files are identical (strict equality after sorting)")
    else:
        print("  ⚠ Strict comparison failed → tolerant numeric check")
        result = compare_numeric_dfs(czid_sorted, seqtoid_sorted, id_cols=[sort_col], wide_matrix=True)
        print("    → " + result)


# Step 5: Per-sample Contig Summary Reports
# (unchanged)

def compare_contig_summary_reports():
    print("\n=== Step 5: Per-sample Contig Summary Reports Comparison ===")
    missing_czid = []
    missing_seqtoid = []

    for sample in EXPECTED_SAMPLES:
        czid_files = glob.glob(os.path.join(CZID_DIR, f"{sample}_*_contig_summary_report.csv"))
        seqtoid_files = glob.glob(os.path.join(SEQTOID_DIR, f"{sample}_*_contig_summary_report.csv"))

        if len(czid_files) != 1:
            missing_czid.append(sample)
            print(f"  ✗ {sample}: missing or multiple contig_summary_report in czid")
            continue
        if len(seqtoid_files) != 1:
            missing_seqtoid.append(sample)
            print(f"  ✗ {sample}: missing or multiple contig_summary_report in seqtoid")
            continue

        czid_path = czid_files[0]
        seqtoid_path = seqtoid_files[0]

        print(f"  → {sample} ...", end=" ")

        try:
            czid_df = pd.read_csv(czid_path, dtype=str)
            seqtoid_df = pd.read_csv(seqtoid_path, dtype=str)

            sort_candidates = ['contig_name', 'contig_id', 'contig', 'name', 'id', 'Contig', 'ContigID']
            sort_col = next((c for c in sort_candidates if c in czid_df.columns), None)

            if sort_col:
                czid_df[sort_col] = czid_df[sort_col].astype(str)
                seqtoid_df[sort_col] = seqtoid_df[sort_col].astype(str)
                czid_sorted = czid_df.sort_values(sort_col).reset_index(drop=True)
                seqtoid_sorted = seqtoid_df.sort_values(sort_col).reset_index(drop=True)
            else:
                czid_sorted = czid_df.reset_index(drop=True)
                seqtoid_sorted = seqtoid_df.reset_index(drop=True)
                print("(no reliable sort column found → comparing in original order)")

            if czid_sorted.equals(seqtoid_sorted):
                print("identical")
            else:
                print("strict fail → tolerant numeric check")
                result = compare_numeric_dfs(czid_sorted, seqtoid_sorted, id_cols=[sort_col] if sort_col else [])
                print(f"      → {result}")

                if len(czid_sorted) != len(seqtoid_sorted):
                    print(f"      Row count: czid = {len(czid_sorted)}, seqtoid = {len(seqtoid_sorted)}")

        except Exception as e:
            print(f"failed to compare: {e}")

    if missing_czid:
        print(f"  Missing in czid: {', '.join(missing_czid)}")
    if missing_seqtoid:
        print(f"  Missing in seqtoid: {', '.join(missing_seqtoid)}")


# ────────────────────────────────────────────────────────────────
# Step 6: Non-host reads comparison (R1 & R2)  ← copied from short-reads version
# ────────────────────────────────────────────────────────────────

def compare_nonhost_fastqs():
    """
    Step 6: Non-host reads comparison (sub-steps 6.1–6.3)
    Compares _reads_nh_R1.fastq and _reads_nh_R2.fastq for each sample.
    """
    print("\n=== Step 6: Non-host reads comparison (R1 & R2) ===")

    for sample in EXPECTED_SAMPLES:
        print(f"\n  Sample: {sample}")

        for read in ['R1', 'R2']:
            pattern = f"{sample}_*_reads_nh_{read}.fastq"

            czid_files = glob.glob(os.path.join(CZID_DIR, pattern))
            seqtoid_files = glob.glob(os.path.join(SEQTOID_DIR, pattern))

            if len(czid_files) != 1 or len(seqtoid_files) != 1:
                status = "missing" if len(czid_files) == 0 else "multiple files"
                print(f"    {read}: {status} in one or both directories")
                continue

            czid_fq = czid_files[0]
            seqtoid_fq = seqtoid_files[0]

            print(f"    {read} FASTQ...")

            # 6.1 – Byte-for-byte identity (fastest check)
            czid_hash = file_sha256(czid_fq)
            seqtoid_hash = file_sha256(seqtoid_fq)

            if czid_hash == seqtoid_hash:
                print("      → Files are byte-for-byte identical")
                continue

            print("      → Files differ byte-for-byte")

            # 6.2 – Same number of reads? (quick line count check)
            czid_lines = sum(1 for _ in open(czid_fq))
            seqtoid_lines = sum(1 for _ in open(seqtoid_fq))

            if czid_lines != seqtoid_lines:
                print(f"      → Different number of lines: czid {czid_lines}, seqtoid {seqtoid_lines}")
                print("        (likely different read count)")
                continue

            # 6.3 – Same sequences (ignore order & qualities) – sorted seq hash
            # Requires seqkit installed and in PATH
            try:
                cmd_czid = f"seqkit seq -s -i {czid_fq} | sort | sha256sum"
                cmd_seqtoid = f"seqkit seq -s -i {seqtoid_fq} | sort | sha256sum"

                czid_seq_hash = subprocess.check_output(cmd_czid, shell=True, text=True).split()[0]
                seqtoid_seq_hash = subprocess.check_output(cmd_seqtoid, shell=True, text=True).split()[0]

                if czid_seq_hash == seqtoid_seq_hash:
                    print("      → Same set of sequences (order & qualities ignored)")
                else:
                    print("      → Different sequences (after sorting)")

            except FileNotFoundError:
                print("      → seqkit not found → skipping sequence hash comparison")
            except subprocess.CalledProcessError as e:
                print(f"      → seqkit failed: {e}")


def main():
    print("CZID Long Reads Pipeline Comparison (czid vs seqtoid)\n")
    compare_metadata()
    compare_sample_overviews()
    compare_taxon_reports()
    compare_combined_taxon_results()
    compare_contig_summary_reports()
    compare_nonhost_fastqs()          # ← new Step 6
    print("\nComparison complete.")


if __name__ == '__main__':
    main()