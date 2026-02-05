import pandas as pd
import os
import glob
import numpy as np
import json
import hashlib
import subprocess
from typing import List, Dict
from scipy.sparse import coo_matrix

# ────────────────────────────────────────────────────────────────
# Configuration
# ────────────────────────────────────────────────────────────────

CZID_DIR = 'czid'
SEQTOID_DIR = 'seqtoid'

EXPECTED_SAMPLES = [
    'ERR11417004',
    'SRR23038836_75M_1',
    'SRR11454628',
    'SRR12876565',
    'SRR10903401',
    'SRR1304850',
    'SRR13227005',
    'SRR13227004',
    'SRR13227003',
    'SRR11278904',
    'SRR14136236_75M'
]

DIFF_SYMBOLS = {
    'equivalent':   '✅ <0.005',
    'warning':      '⚠️ 0.005–0.05',
    'significant':  '❌ >0.05',
    'no diffs':     'identical',
    'identical':    'T',
    'differ':       'F'
}

# ────────────────────────────────────────────────────────────────
# Helpers
# ────────────────────────────────────────────────────────────────

def file_sha256(filepath):
    """Compute SHA-256 hash of a file in chunks."""
    sha256 = hashlib.sha256()
    with open(filepath, 'rb') as f:
        while chunk := f.read(8192 * 1024):
            sha256.update(chunk)
    return sha256.hexdigest()


def numeric_diff(a: np.ndarray, b: np.ndarray, atol: float = 0.005) -> str:
    """Categorize worst difference."""
    if len(a) == 0 or len(b) == 0:
        return "empty"
    diff = np.abs(a - b)
    worst = diff.max()
    if np.isnan(worst):
        return "NaN"
    if worst <= 0.005:
        return 'equivalent'
    elif worst <= 0.05:
        return 'warning'
    else:
        return 'significant'


def compare_numeric_dfs(df1: pd.DataFrame, df2: pd.DataFrame, id_cols: List[str]) -> Dict[str, str]:
    """Return dict of numeric column → category."""
    num_cols = df1.select_dtypes(include=[np.number]).columns.intersection(df2.columns)
    results = {}
    for col in num_cols:
        vals1 = df1[col].fillna(0).values.astype(float)
        vals2 = df2[col].fillna(0).values.astype(float)
        cat = numeric_diff(vals1, vals2)
        results[col] = cat
    return results


# ────────────────────────────────────────────────────────────────
# Step functions with CSV output
# ────────────────────────────────────────────────────────────────

def compare_metadata():
    file = 'sample_metadata.csv'
    czid_path = os.path.join(CZID_DIR, file)
    seqtoid_path = os.path.join(SEQTOID_DIR, file)

    if not os.path.exists(czid_path) or not os.path.exists(seqtoid_path):
        print(f"Error: {file} missing in one or both dirs.")
        pd.DataFrame({'file': [file], 'identical': ['missing']}).to_csv('step1_metadata.csv', index=False)
        return

    czid_df = pd.read_csv(czid_path)
    seqtoid_df = pd.read_csv(seqtoid_path)

    czid_s = czid_df.sort_values('sample_name').reset_index(drop=True)
    seqtoid_s = seqtoid_df.sort_values('sample_name').reset_index(drop=True)

    identical = czid_s.equals(seqtoid_s)
    print(f"Step 1: sample_metadata.csv → {'identical' if identical else 'differ'}")

    pd.DataFrame({'file': [file], 'identical': ['T' if identical else 'F']}).to_csv('step1_metadata.csv', index=False)

    # Missing / extra
    rows = []
    for name, df in [('czid', czid_df), ('seqtoid', seqtoid_df)]:
        actual = set(df['sample_name'].astype(str))
        missing = ', '.join(sorted(set(EXPECTED_SAMPLES) - actual))
        extra = ', '.join(sorted(actual - set(EXPECTED_SAMPLES)))
        if missing or extra:
            rows.append({'directory': name, 'missing': missing, 'extra': extra})
    if rows:
        pd.DataFrame(rows).to_csv('step1_missing_extra.csv', index=False)


def compare_overviews():
    file_name = 'sample_overviews.csv'
    czid_path = os.path.join(CZID_DIR, file_name)
    seqtoid_path = os.path.join(SEQTOID_DIR, file_name)

    if not os.path.exists(czid_path) or not os.path.exists(seqtoid_path):
        print(f"Error: One or both {file_name} files missing.")
        pd.DataFrame({'file': [file_name], 'status': ['missing']}).to_csv('step2_overviews_comparison.csv', index=False)
        return

    czid_df = pd.read_csv(czid_path, dtype=str)
    seqtoid_df = pd.read_csv(seqtoid_path, dtype=str)

    czid_sorted = czid_df.sort_values('sample_name').reset_index(drop=True)
    seqtoid_sorted = seqtoid_df.sort_values('sample_name').reset_index(drop=True)

    print("Comparing sample_overviews.csv...")

    # Always perform numeric comparison
    print("  → Performing numeric tolerance check...")

    # Convert candidate numeric columns to float
    numeric_candidates = [
        'runtime_seconds', 'total_reads', 'passed_filters', 'passed_filters_percent',
        'total_ercc_reads', 'subsampled_fraction', 'quality_control', 'compression_ratio',
        'reads_after_bowtie2_ercc_filtered', 'reads_after_fastp',
        'reads_after_bowtie2_host_filtered', 'reads_after_hisat2_host_filtered',
        'reads_after_czid_dedup', 'insert_size_median', 'insert_size_mode',
        'insert_size_median_absolute_deviation', 'insert_size_min', 'insert_size_max',
        'insert_size_mean', 'insert_size_standard_deviation', 'insert_size_read_pairs'
    ]

    for col in numeric_candidates:
        if col in czid_sorted.columns:
            czid_sorted[col] = pd.to_numeric(czid_sorted[col], errors='coerce')
            seqtoid_sorted[col] = pd.to_numeric(seqtoid_sorted[col], errors='coerce')

    # Compare only numeric columns
    diffs = compare_numeric_dfs(czid_sorted, seqtoid_sorted, id_cols=['sample_name'])

    # Prepare results
    rows = []
    if diffs:
        for col, cat in diffs.items():
            rows.append({
                'column': col,
                'category': cat,
                'symbol': DIFF_SYMBOLS.get(cat, cat)
            })
    else:
        rows.append({'status': 'no numeric differences detected'})

    pd.DataFrame(rows).to_csv('step2_overviews_comparison.csv', index=False)
    print("  → Results written to step2_overviews_comparison.csv")


def compare_taxon_reports():
    print("Step 3: per-sample taxon_report.csv")
    rows = []
    missing_czid = []
    missing_seqtoid = []

    for sample in EXPECTED_SAMPLES:
        pattern = f"{sample}_*_taxon_report.csv"
        czid_f = glob.glob(os.path.join(CZID_DIR, pattern))
        seqtoid_f = glob.glob(os.path.join(SEQTOID_DIR, pattern))

        if len(czid_f) != 1 or len(seqtoid_f) != 1:
            if len(czid_f) != 1: missing_czid.append(sample)
            if len(seqtoid_f) != 1: missing_seqtoid.append(sample)
            rows.append({'sample': sample, 'status': 'missing'})
            continue

        czid_df = pd.read_csv(czid_f[0], dtype=str)
        seqtoid_df = pd.read_csv(seqtoid_f[0], dtype=str)

        czid_df['tax_id'] = pd.to_numeric(czid_df['tax_id'], errors='coerce')
        seqtoid_df['tax_id'] = pd.to_numeric(seqtoid_df['tax_id'], errors='coerce')

        czid_s = czid_df.sort_values('tax_id').reset_index(drop=True)
        seqtoid_s = seqtoid_df.sort_values('tax_id').reset_index(drop=True)

        if czid_s.equals(seqtoid_s):
            rows.append({'sample': sample, 'status': 'identical'})
        else:
            diffs = compare_numeric_dfs(czid_s, seqtoid_s, id_cols=['tax_id'])
            for col, cat in diffs.items():
                rows.append({'sample': sample, 'column': col, 'category': cat, 'symbol': DIFF_SYMBOLS.get(cat, cat)})

    pd.DataFrame(rows).to_csv('step3_taxon_reports.csv', index=False)

    if missing_czid or missing_seqtoid:
        pd.DataFrame({
            'missing_in_czid': [', '.join(missing_czid)],
            'missing_in_seqtoid': [', '.join(missing_seqtoid)]
        }).to_csv('step3_missing.csv', index=False)


def compare_combined_taxon_results():
    """
    Step 4: Compare combined_sample_taxon_results_NT.rpm.csv
    - Aligns on 'Taxon Name'
    - Per-taxon per-sample cell comparison:
      - Both values: numeric diff category
      - One value, one missing/NaN: 'missing'
      - Both missing/NaN: skip
    - Writes long-form CSV with per-cell results
    """
    file_name = 'combined_sample_taxon_results_NT.rpm.csv'
    czid_path = os.path.join(CZID_DIR, file_name)
    seqtoid_path = os.path.join(SEQTOID_DIR, file_name)

    if not os.path.exists(czid_path) or not os.path.exists(seqtoid_path):
        print(f"Error: One or both {file_name} files missing.")
        pd.DataFrame({'file': [file_name], 'status': ['missing']}).to_csv('step4_combined_rpm.csv', index=False)
        return

    czid_df = pd.read_csv(czid_path)
    seqtoid_df = pd.read_csv(seqtoid_path)

    # Assume first column is 'Taxon Name'
    taxon_col = czid_df.columns[0]

    # Set index to taxon for merging
    czid_df.set_index(taxon_col, inplace=True)
    seqtoid_df.set_index(taxon_col, inplace=True)

    # Outer merge to include all taxa
    merged = czid_df.join(seqtoid_df, lsuffix='_czid', rsuffix='_seqtoid', how='outer')

    # Sample columns (assume all except index are samples)
    sample_cols = [col.replace('_czid', '') for col in merged.columns if '_czid' in col]

    # Collect per-cell results
    rows = []
    for taxon in merged.index:
        for sample in sample_cols:
            val_czid = merged.at[taxon, f"{sample}_czid"]
            val_seqtoid = merged.at[taxon, f"{sample}_seqtoid"]

            # Skip if both NaN/missing
            if pd.isna(val_czid) and pd.isna(val_seqtoid):
                continue

            # Missing if one NaN and other not
            if pd.isna(val_czid) or pd.isna(val_seqtoid):
                category = 'missing'
                abs_diff = 'N/A'
            else:
                # Numeric comparison
                val_czid = float(val_czid)
                val_seqtoid = float(val_seqtoid)
                abs_diff = abs(val_czid - val_seqtoid)
                if abs_diff <= 0.005:
                    category = 'equivalent'
                elif abs_diff <= 0.05:
                    category = 'warning'
                else:
                    category = 'significant'

            symbol = DIFF_SYMBOLS.get(category, category)

            rows.append({
                'taxon': taxon,
                'sample': sample,
                'value_czid': val_czid,
                'value_seqtoid': val_seqtoid,
                'abs_diff': abs_diff if category != 'missing' else 'N/A',
                'category': category,
                'symbol': symbol
            })

    # Write long-form CSV
    pd.DataFrame(rows).to_csv('step4_combined_rpm.csv', index=False)
    print(f"  → Per-cell comparison complete. Results in step4_combined_rpm.csv ({len(rows)} entries)")


def compare_contig_summary_reports():
    """
    Step 5: Compare per-sample contig_summary_report.csv files.
    - For each sample, matches contig_name rows and checks if NT.species_taxid matches.
    - Counts missing contigs as mismatch.
    - Computes mismatch proportion (non-matching NT.species_taxid / total unique contig_names).
    - Categories: <0.05 equivalent, 0.05-0.1 warning, >0.1 significant.
    - Writes per-sample results to CSV.
    """
    print("Comparing per-sample contig_summary_report.csv files...")

    results = []
    missing_in_czid = []
    missing_in_seqtoid = []

    for sample in EXPECTED_SAMPLES:
        # Find files
        czid_pattern = os.path.join(CZID_DIR, f"{sample}_*_contig_summary_report.csv")
        seqtoid_pattern = os.path.join(SEQTOID_DIR, f"{sample}_*_contig_summary_report.csv")

        czid_files = glob.glob(czid_pattern)
        seqtoid_files = glob.glob(seqtoid_pattern)

        if len(czid_files) != 1 or len(seqtoid_files) != 1:
            print(f"  - For sample {sample}: {'No file found' if len(czid_files) == 0 else 'Multiple files'} in czid")
            print(f"  - For sample {sample}: {'No file found' if len(seqtoid_files) == 0 else 'Multiple files'} in seqtoid")
            if len(czid_files) != 1:
                missing_in_czid.append(sample)
            if len(seqtoid_files) != 1:
                missing_in_seqtoid.append(sample)
            results.append({'sample': sample, 'mismatch_proportion': 'N/A', 'category': 'missing', 'symbol': 'missing'})
            continue

        czid_path = czid_files[0]
        seqtoid_path = seqtoid_files[0]

        # Read CSVs
        czid_df = pd.read_csv(czid_path)
        seqtoid_df = pd.read_csv(seqtoid_path)

        # Use contig_name as index for easy lookup
        czid_df.set_index('contig_name', inplace=True)
        seqtoid_df.set_index('contig_name', inplace=True)

        # Get all unique contig_names from both
        all_contigs = set(czid_df.index) | set(seqtoid_df.index)
        total_contigs = len(all_contigs)

        # Count mismatches on NT.species_taxid
        mismatches = 0
        mismatched_list = []  # Optional: collect mismatched contigs

        for contig in all_contigs:
            if contig not in czid_df.index or contig not in seqtoid_df.index:
                mismatches += 1
                mismatched_list.append(contig)
                continue

            # Compare NT.species_taxid (convert to string to handle NaN/-100 etc. safely)
            czid_taxid = str(czid_df.at[contig, 'NT.species_taxid']) if 'NT.species_taxid' in czid_df.columns else None
            seqtoid_taxid = str(seqtoid_df.at[contig, 'NT.species_taxid']) if 'NT.species_taxid' in seqtoid_df.columns else None

            if czid_taxid != seqtoid_taxid:
                mismatches += 1
                mismatched_list.append(contig)

        # Compute proportion
        proportion = mismatches / total_contigs if total_contigs > 0 else 0

        # Categorize
        if proportion < 0.05:
            category = 'equivalent'
        elif proportion < 0.1:
            category = 'warning'
        else:
            category = 'significant'

        symbol = DIFF_SYMBOLS.get(category, category)

        print(f"  - {sample}: mismatch proportion {proportion:.4f} ({mismatches}/{total_contigs}) → {symbol}")

        results.append({
            'sample': sample,
            'total_contigs': total_contigs,
            'mismatches': mismatches,
            'mismatch_proportion': round(proportion, 6),
            'category': category,
            'symbol': symbol
        })

    # Write overall CSV
    pd.DataFrame(results).to_csv('step5_contig_summary_comparison.csv', index=False)

    # Optional: if you want to see which contigs mismatched for debugging
    # You could write mismatched_list to a separate file per sample, but keeping it simple here

    if missing_in_czid:
        print(f"  - Missing in czid: {', '.join(missing_in_czid)}")
    if missing_in_seqtoid:
        print(f"  - Missing in seqtoid: {', '.join(missing_in_seqtoid)}")


def compare_host_gene_counts():
    """
    Step 6: Compare per-sample reads_per_transcript.kallisto.tsv files.
    - For each sample, counts unique target_id in czid and seqtoid.
    - Computes common, unique to czid, unique to seqtoid.
    - Proportion discrepancies = (unique_czid + unique_seqtoid) / total_unique.
    - Categories based on proportion: <0.05 equivalent, 0.05-0.1 warning, >0.1 significant.
    - Writes per-sample results to CSV.
    """
    print("Comparing per-sample reads_per_transcript.kallisto.tsv files...")

    results = []
    missing_in_czid = []
    missing_in_seqtoid = []

    for sample in EXPECTED_SAMPLES:
        # Find files
        czid_pattern = os.path.join(CZID_DIR, f"{sample}_*_reads_per_transcript.kallisto.tsv")
        seqtoid_pattern = os.path.join(SEQTOID_DIR, f"{sample}_*_reads_per_transcript.kallisto.tsv")

        czid_files = glob.glob(czid_pattern)
        seqtoid_files = glob.glob(seqtoid_pattern)

        if len(czid_files) != 1 or len(seqtoid_files) != 1:
            print(f"  - For sample {sample}: {'No file found' if len(czid_files) == 0 else 'Multiple files'} in czid")
            print(f"  - For sample {sample}: {'No file found' if len(seqtoid_files) == 0 else 'Multiple files'} in seqtoid")
            if len(czid_files) != 1:
                missing_in_czid.append(sample)
            if len(seqtoid_files) != 1:
                missing_in_seqtoid.append(sample)
            results.append({'sample': sample, 'discrepancy_proportion': 'N/A', 'category': 'missing', 'symbol': 'missing'})
            continue

        czid_path = czid_files[0]
        seqtoid_path = seqtoid_files[0]

        # Read TSVs
        czid_df = pd.read_csv(czid_path, sep='\t')
        seqtoid_df = pd.read_csv(seqtoid_path, sep='\t')

        # Get unique target_id sets
        czid_ids = set(czid_df['target_id'].unique())
        seqtoid_ids = set(seqtoid_df['target_id'].unique())

        # Counts
        unique_czid = len(czid_ids - seqtoid_ids)
        unique_seqtoid = len(seqtoid_ids - czid_ids)
        common = len(czid_ids & seqtoid_ids)
        total_unique = unique_czid + unique_seqtoid + common

        # Proportion of discrepancies (unique to either / total unique)
        discrepancy_proportion = (unique_czid + unique_seqtoid) / total_unique if total_unique > 0 else 0

        # Categorize
        if discrepancy_proportion < 0.05:
            category = 'equivalent'
        elif discrepancy_proportion < 0.1:
            category = 'warning'
        else:
            category = 'significant'

        symbol = DIFF_SYMBOLS.get(category, category)

        print(f"  - {sample}: discrepancy proportion {discrepancy_proportion:.4f} (unique czid {unique_czid}, unique seqtoid {unique_seqtoid}, common {common}) → {symbol}")

        results.append({
            'sample': sample,
            'unique_czid': unique_czid,
            'unique_seqtoid': unique_seqtoid,
            'common': common,
            'total_unique': total_unique,
            'discrepancy_proportion': round(discrepancy_proportion, 6),
            'category': category,
            'symbol': symbol
        })

    # Write CSV
    pd.DataFrame(results).to_csv('step6_host_counts_comparison.csv', index=False)

    if missing_in_czid:
        print(f"  - Missing in czid: {', '.join(missing_in_czid)}")
    if missing_in_seqtoid:
        print(f"  - Missing in seqtoid: {', '.join(missing_in_seqtoid)}")


def compare_combined_microbiome():
    file = 'Combined Microbiome File.biom'
    czid_path = os.path.join(CZID_DIR, file)
    seqtoid_path = os.path.join(SEQTOID_DIR, file)

    if not os.path.exists(czid_path) or not os.path.exists(seqtoid_path):
        print(f"Error: {file} missing in one or both dirs.")
        pd.DataFrame({'file': [file], 'status': ['missing']}).to_csv('step7_biom.csv', index=False)
        return

    # Load both
    with open(czid_path, 'r') as f:
        czid = json.load(f)
    with open(seqtoid_path, 'r') as f:
        seqtoid = json.load(f)

    print("Step 7: BIOM file")

    # Very simplified for CSV – extend as needed
    rows = [{'file': file, 'status': 'compared'}]
    # Add more rows if you expand the comparison (e.g. shape match, row count match, value diff category)

    pd.DataFrame(rows).to_csv('step7_biom.csv', index=False)


def compare_nonhost_fastqs():
    print("Step 8: Non-host reads FASTQ (R1 & R2)")
    rows = []
    missing_czid = []
    missing_seqtoid = []

    for sample in EXPECTED_SAMPLES:
        for read in ['R1', 'R2']:
            pattern = f"{sample}_*_reads_nh_{read}.fastq"
            czid_f = glob.glob(os.path.join(CZID_DIR, pattern))
            seqtoid_f = glob.glob(os.path.join(SEQTOID_DIR, pattern))

            if len(czid_f) != 1 or len(seqtoid_f) != 1:
                if len(czid_f) != 1: missing_czid.append(f"{sample} {read}")
                if len(seqtoid_f) != 1: missing_seqtoid.append(f"{sample} {read}")
                rows.append({'sample': sample, 'read': read, 'identical': 'missing'})
                continue

            czid_fa = czid_f[0]
            seqtoid_fa = seqtoid_f[0]

            identical = file_sha256(czid_fa) == file_sha256(seqtoid_fa)
            rows.append({'sample': sample, 'read': read, 'identical': 'T' if identical else 'F'})

    pd.DataFrame(rows).to_csv('step8_nonhost_fastq.csv', index=False)

    if missing_czid or missing_seqtoid:
        pd.DataFrame({
            'missing_czid': [', '.join(missing_czid)],
            'missing_seqtoid': [', '.join(missing_seqtoid)]
        }).to_csv('step8_missing.csv', index=False)


def compare_nonhost_contigs():
    print("Step 9: Non-host contigs FASTA")
    rows = []
    missing_czid = []
    missing_seqtoid = []

    for sample in EXPECTED_SAMPLES:
        pattern = f"{sample}_*_contigs_nh.fasta"
        czid_f = glob.glob(os.path.join(CZID_DIR, pattern))
        seqtoid_f = glob.glob(os.path.join(SEQTOID_DIR, pattern))

        if len(czid_f) != 1 or len(seqtoid_f) != 1:
            if len(czid_f) != 1: missing_czid.append(sample)
            if len(seqtoid_f) != 1: missing_seqtoid.append(sample)
            rows.append({'sample': sample, 'identical': 'missing'})
            continue

        czid_fa = czid_f[0]
        seqtoid_fa = seqtoid_f[0]

        identical = file_sha256(czid_fa) == file_sha256(seqtoid_fa)
        rows.append({'sample': sample, 'identical': 'T' if identical else 'F'})

    pd.DataFrame(rows).to_csv('step9_nonhost_contigs.csv', index=False)

    if missing_czid or missing_seqtoid:
        pd.DataFrame({
            'missing_in_czid': [', '.join(missing_czid)],
            'missing_in_seqtoid': [', '.join(missing_seqtoid)]
        }).to_csv('step9_missing.csv', index=False)


def main():
    print("Starting CZID short-read pipeline comparison...\n")

    print("=== Step 1: Sample Metadata Comparison ===")
    compare_metadata()

    print("\n=== Step 2: Sample Overviews Comparison ===")
    compare_overviews()

    print("\n=== Step 3: Sample Taxon Reports Comparison ===")
    compare_taxon_reports()

    print("\n=== Step 4: Combined Taxon RPM Comparison ===")
    compare_combined_taxon_results()

    print("\n=== Step 5: Contig Summary Reports Comparison ===")
    compare_contig_summary_reports()

    print("\n=== Step 6: Host Gene Counts (Kallisto) Comparison ===")
    compare_host_gene_counts()

    print("\n=== Step 7: Combined Microbiome BIOM Comparison ===")
    compare_combined_microbiome()

    print("\n=== Step 8: Non-host reads FASTQ Comparison ===")
    compare_nonhost_fastqs()

    print("\n=== Step 9: Non-host contigs FASTA Comparison ===")
    compare_nonhost_contigs()

    print("\nComparison complete. Check CSV files in current directory.")


if __name__ == '__main__':
    main()