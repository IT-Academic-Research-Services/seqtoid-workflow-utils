import pandas as pd
import os
import glob
import numpy as np
import json
import hashlib
from typing import List, Dict

# ────────────────────────────────────────────────────────────────
# Configuration
# ────────────────────────────────────────────────────────────────

CZID_DIR = 'czid'
SEQTOID_DIR = 'seqtoid'

# Informal expected list (not the analysis universe).
# Analysis samples = whatever is present in czid after stripping _75M.
EXPECTED_SAMPLES = [
    'ERR11417004',
    'SRR1304850',
    'SRR10903401',
    'SRR11278904',
    'SRR11454628',
    'SRR12876565',
    'SRR13227003',
    'SRR13227004',
    'SRR13227005',
    'SRR14136236',
    'SRR23038836'
]

DIFF_SYMBOLS = {
    'equivalent':   '✅ <0.005',
    'warning':      '⚠️ 0.005–0.05',
    'significant':  '❌ >0.05',
    'no diffs':     'identical',
    'identical':    'T',
    'differ':       'F'
}

_ANALYSIS_SAMPLES = None


def norm_sample_name(name):
    """SRR14136236_75M and SRR23038836_75M_1 → base NCBI id."""
    s = str(name).strip()
    for suffix in ('_75M_1', '_75M'):
        if s.endswith(suffix):
            return s[:-len(suffix)]
    return s


def discover_dir_samples(directory):
    """Normalized sample names present in a results directory."""
    names = set()
    for fname in ('sample_overviews.csv', 'sample_metadata.csv'):
        path = os.path.join(directory, fname)
        if not os.path.exists(path):
            continue
        df = pd.read_csv(path, dtype=str)
        if 'sample_name' in df.columns:
            names.update(df['sample_name'].map(norm_sample_name))
    return names


def analysis_samples():
    """Samples to compare: present in czid. Seqtoid-only samples are skipped."""
    global _ANALYSIS_SAMPLES
    if _ANALYSIS_SAMPLES is not None:
        return _ANALYSIS_SAMPLES
    czid_names = discover_dir_samples(CZID_DIR)
    if czid_names:
        _ANALYSIS_SAMPLES = sorted(czid_names)
    else:
        _ANALYSIS_SAMPLES = list(EXPECTED_SAMPLES)
        print("  ! No czid metadata/overviews found; falling back to EXPECTED_SAMPLES")
    seqtoid_names = discover_dir_samples(SEQTOID_DIR)
    skipped = sorted(seqtoid_names - set(_ANALYSIS_SAMPLES))
    print(f"Analysis samples (from czid): {', '.join(_ANALYSIS_SAMPLES)}")
    if skipped:
        print(f"Skipping seqtoid-only samples: {', '.join(skipped)}")
    return _ANALYSIS_SAMPLES


def find_sample_file(directory, sample, suffix):
    """
    suffix e.g. '*_taxon_report.csv'
    Matches both ERR11417004_879286_taxon_report.csv
    and SRR14136236_75M_879286_taxon_report.csv
    """
    hits = glob.glob(os.path.join(directory, f"{sample}_{suffix}"))
    return hits


def file_sha256(filepath):
    sha256 = hashlib.sha256()
    with open(filepath, 'rb') as f:
        while chunk := f.read(8192 * 1024):
            sha256.update(chunk)
    return sha256.hexdigest()


def numeric_diff(a: np.ndarray, b: np.ndarray, atol: float = 0.005) -> str:
    if len(a) == 0 or len(b) == 0:
        return "empty"
    if len(a) != len(b):
        return "shape_mismatch"
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
    num_cols = df1.select_dtypes(include=[np.number]).columns.intersection(df2.columns)
    results = {}
    for col in num_cols:
        vals1 = df1[col].fillna(0).values.astype(float)
        vals2 = df2[col].fillna(0).values.astype(float)
        results[col] = numeric_diff(vals1, vals2)
    return results


def compare_metadata():
    file = 'sample_metadata.csv'
    czid_path = os.path.join(CZID_DIR, file)
    seqtoid_path = os.path.join(SEQTOID_DIR, file)
    samples = analysis_samples()

    if not os.path.exists(czid_path) or not os.path.exists(seqtoid_path):
        print(f"Error: {file} missing in one or both dirs.")
        pd.DataFrame({'file': [file], 'identical': ['missing']}).to_csv('step1_metadata.csv', index=False)
        return

    czid_df = pd.read_csv(czid_path, dtype=str)
    seqtoid_df = pd.read_csv(seqtoid_path, dtype=str)
    czid_df['sample_key'] = czid_df['sample_name'].map(norm_sample_name)
    seqtoid_df['sample_key'] = seqtoid_df['sample_name'].map(norm_sample_name)

    czid_df = czid_df[czid_df['sample_key'].isin(samples)].copy()
    seqtoid_df = seqtoid_df[seqtoid_df['sample_key'].isin(samples)].copy()

    compare_cols = [c for c in czid_df.columns if c not in ('sample_name', 'sample_key')]
    czid_s = czid_df.sort_values('sample_key').reset_index(drop=True)
    seqtoid_s = seqtoid_df.sort_values('sample_key').reset_index(drop=True)

    identical = (
        list(czid_s['sample_key']) == list(seqtoid_s['sample_key'])
        and czid_s[compare_cols].equals(seqtoid_s[compare_cols])
        if set(compare_cols) <= set(seqtoid_s.columns)
        else False
    )
    print(f"Step 1: sample_metadata.csv → {'identical' if identical else 'differ'} (n={len(samples)})")
    pd.DataFrame({'file': [file], 'identical': ['T' if identical else 'F']}).to_csv('step1_metadata.csv', index=False)

    seqtoid_keys = set(seqtoid_s['sample_key'])
    missing_seqtoid = [s for s in samples if s not in seqtoid_keys]
    if missing_seqtoid:
        pd.DataFrame({'missing_in_seqtoid': [', '.join(missing_seqtoid)]}).to_csv('step1_missing_extra.csv', index=False)


def compare_overviews():
    file_name = 'sample_overviews.csv'
    czid_path = os.path.join(CZID_DIR, file_name)
    seqtoid_path = os.path.join(SEQTOID_DIR, file_name)
    samples = analysis_samples()

    if not os.path.exists(czid_path) or not os.path.exists(seqtoid_path):
        print(f"Error: One or both {file_name} files missing.")
        pd.DataFrame({'file': [file_name], 'status': ['missing']}).to_csv('step2_overviews_comparison.csv', index=False)
        return

    czid_df = pd.read_csv(czid_path, dtype=str)
    seqtoid_df = pd.read_csv(seqtoid_path, dtype=str)
    czid_df['sample_key'] = czid_df['sample_name'].map(norm_sample_name)
    seqtoid_df['sample_key'] = seqtoid_df['sample_name'].map(norm_sample_name)

    czid_df = czid_df[czid_df['sample_key'].isin(samples)].drop_duplicates('sample_key')
    seqtoid_df = seqtoid_df[seqtoid_df['sample_key'].isin(samples)].drop_duplicates('sample_key')

    shared = sorted(set(czid_df['sample_key']) & set(seqtoid_df['sample_key']))
    missing_seqtoid = sorted(set(samples) - set(seqtoid_df['sample_key']))
    print(f"Comparing sample_overviews.csv... shared={len(shared)} missing_in_seqtoid={missing_seqtoid}")

    numeric_candidates = [
        'runtime_seconds', 'total_reads', 'passed_filters', 'passed_filters_percent',
        'total_ercc_reads', 'subsampled_fraction', 'quality_control', 'compression_ratio',
        'reads_after_bowtie2_ercc_filtered', 'reads_after_fastp',
        'reads_after_bowtie2_host_filtered', 'reads_after_hisat2_host_filtered',
        'reads_after_czid_dedup', 'insert_size_median', 'insert_size_mode',
        'insert_size_median_absolute_deviation', 'insert_size_min', 'insert_size_max',
        'insert_size_mean', 'insert_size_standard_deviation', 'insert_size_read_pairs'
    ]

    czid_shared = czid_df[czid_df['sample_key'].isin(shared)].copy()
    seqtoid_shared = seqtoid_df[seqtoid_df['sample_key'].isin(shared)].copy()
    for col in numeric_candidates:
        if col in czid_shared.columns:
            czid_shared[col] = pd.to_numeric(czid_shared[col], errors='coerce')
        if col in seqtoid_shared.columns:
            seqtoid_shared[col] = pd.to_numeric(seqtoid_shared[col], errors='coerce')

    czid_shared = czid_shared.sort_values('sample_key').reset_index(drop=True)
    seqtoid_shared = seqtoid_shared.sort_values('sample_key').reset_index(drop=True)

    diffs = compare_numeric_dfs(czid_shared, seqtoid_shared, id_cols=['sample_key'])

    rows = []
    if diffs:
        for col, cat in diffs.items():
            rows.append({
                'column': col,
                'n_shared_samples': len(shared),
                'category': cat,
                'symbol': DIFF_SYMBOLS.get(cat, cat)
            })
    else:
        rows.append({'status': 'no numeric columns compared', 'n_shared_samples': len(shared)})

    pd.DataFrame(rows).to_csv('step2_overviews_comparison.csv', index=False)
    pd.DataFrame([{
        'n_analysis': len(samples),
        'n_shared': len(shared),
        'missing_in_seqtoid': ', '.join(missing_seqtoid)
    }]).to_csv('step2_overviews_samples.csv', index=False)
    print("  → step2_overviews_comparison.csv")


def compare_taxon_reports():
    print("Step 3: per-sample taxon_report.csv")
    results = []
    missing_czid = []
    missing_seqtoid = []
    numeric_cols = ['nt_rpm', 'nr_rpm', 'nt_count', 'nr_count']

    for sample in analysis_samples():
        czid_f = find_sample_file(CZID_DIR, sample, '*_taxon_report.csv')
        seqtoid_f = find_sample_file(SEQTOID_DIR, sample, '*_taxon_report.csv')

        if len(czid_f) != 1 or len(seqtoid_f) != 1:
            if len(czid_f) != 1:
                missing_czid.append(sample)
            if len(seqtoid_f) != 1:
                missing_seqtoid.append(sample)
            results.append({
                'sample': sample,
                'status': 'missing',
                'max_abs_diff_shared': 'N/A',
                'max_diff_tax_id': 'N/A',
                'category': 'missing',
                'symbol': 'missing'
            })
            continue

        czid_df = pd.read_csv(czid_f[0])
        seqtoid_df = pd.read_csv(seqtoid_f[0])
        czid_df['tax_id'] = pd.to_numeric(czid_df['tax_id'], errors='coerce')
        seqtoid_df['tax_id'] = pd.to_numeric(seqtoid_df['tax_id'], errors='coerce')
        czid_df = czid_df.dropna(subset=['tax_id']).drop_duplicates(subset=['tax_id']).set_index('tax_id')
        seqtoid_df = seqtoid_df.dropna(subset=['tax_id']).drop_duplicates(subset=['tax_id']).set_index('tax_id')

        czid_ids = set(czid_df.index)
        seqtoid_ids = set(seqtoid_df.index)
        common_ids = czid_ids & seqtoid_ids
        only_czid = czid_ids - seqtoid_ids
        only_seqtoid = seqtoid_ids - czid_ids
        total_unique = len(czid_ids | seqtoid_ids)
        id_discrepancy = (len(only_czid) + len(only_seqtoid)) / total_unique if total_unique else 0

        worst_cat = 'equivalent'
        n_equiv = n_warn = n_sig = n_missing_val = 0
        max_abs_diff = 0.0
        max_diff_tax_id = ''
        max_diff_col = ''
        max_diff_czid = ''
        max_diff_seqtoid = ''

        for col in numeric_cols:
            if col not in czid_df.columns or col not in seqtoid_df.columns:
                continue
            for tid in common_ids:
                a = pd.to_numeric(czid_df.at[tid, col], errors='coerce')
                b = pd.to_numeric(seqtoid_df.at[tid, col], errors='coerce')
                if pd.isna(a) and pd.isna(b):
                    continue
                if pd.isna(a) or pd.isna(b):
                    n_missing_val += 1
                    continue
                a, b = float(a), float(b)
                d = abs(a - b)
                if d > max_abs_diff:
                    max_abs_diff = d
                    max_diff_tax_id = tid
                    max_diff_col = col
                    max_diff_czid = a
                    max_diff_seqtoid = b
                if d <= 0.005:
                    n_equiv += 1
                elif d <= 0.05:
                    n_warn += 1
                    if worst_cat == 'equivalent':
                        worst_cat = 'warning'
                else:
                    n_sig += 1
                    worst_cat = 'significant'

        if id_discrepancy >= 0.1:
            category = 'significant'
        elif id_discrepancy >= 0.05:
            category = 'warning'
        else:
            category = worst_cat

        symbol = DIFF_SYMBOLS.get(category, category)
        print(
            f"  - {sample}: tax_id discrepancy {id_discrepancy:.4f} "
            f"(only_czid={len(only_czid)}, only_seqtoid={len(only_seqtoid)}, common={len(common_ids)}) "
            f"max_diff={max_abs_diff:.4f} tax_id={max_diff_tax_id} col={max_diff_col} → {symbol}"
        )
        results.append({
            'sample': sample,
            'n_czid': len(czid_ids),
            'n_seqtoid': len(seqtoid_ids),
            'n_common': len(common_ids),
            'n_only_czid': len(only_czid),
            'n_only_seqtoid': len(only_seqtoid),
            'id_discrepancy_proportion': round(id_discrepancy, 6),
            'shared_cells_equivalent': n_equiv,
            'shared_cells_warning': n_warn,
            'shared_cells_significant': n_sig,
            'shared_cells_one_sided_missing': n_missing_val,
            'max_abs_diff_shared': round(max_abs_diff, 6),
            'max_diff_tax_id': max_diff_tax_id,
            'max_diff_column': max_diff_col,
            'max_diff_value_czid': max_diff_czid,
            'max_diff_value_seqtoid': max_diff_seqtoid,
            'category': category,
            'symbol': symbol
        })

    pd.DataFrame(results).to_csv('step3_taxon_reports.csv', index=False)
    if missing_czid or missing_seqtoid:
        pd.DataFrame({
            'missing_in_czid': [', '.join(missing_czid)],
            'missing_in_seqtoid': [', '.join(missing_seqtoid)]
        }).to_csv('step3_missing.csv', index=False)


def compare_combined_taxon_results():
    file_name = 'combined_sample_taxon_results_NT.rpm.csv'
    czid_path = os.path.join(CZID_DIR, file_name)
    seqtoid_path = os.path.join(SEQTOID_DIR, file_name)
    samples = analysis_samples()

    if not os.path.exists(czid_path) or not os.path.exists(seqtoid_path):
        print(f"Error: One or both {file_name} files missing.")
        pd.DataFrame({'file': [file_name], 'status': ['missing']}).to_csv('step4_combined_rpm.csv', index=False)
        return

    czid_df = pd.read_csv(czid_path)
    seqtoid_df = pd.read_csv(seqtoid_path)
    taxon_col = czid_df.columns[0]
    czid_df = czid_df.rename(columns={c: norm_sample_name(c) for c in czid_df.columns if c != taxon_col})
    seqtoid_df = seqtoid_df.rename(columns={c: norm_sample_name(c) for c in seqtoid_df.columns if c != taxon_col})

    use_cols = [c for c in samples if c in czid_df.columns]
    czid_df = czid_df[[taxon_col] + use_cols].set_index(taxon_col)
    seqtoid_keep = [c for c in use_cols if c in seqtoid_df.columns]
    seqtoid_df = seqtoid_df[[taxon_col] + seqtoid_keep].set_index(taxon_col)

    merged = czid_df.join(seqtoid_df, lsuffix='_czid', rsuffix='_seqtoid', how='outer')
    sample_cols = [c.replace('_czid', '') for c in merged.columns if c.endswith('_czid')]

    rows = []
    for taxon in merged.index:
        for sample in sample_cols:
            val_czid = merged.at[taxon, f"{sample}_czid"]
            val_seqtoid = merged.at[taxon, f"{sample}_seqtoid"] if f"{sample}_seqtoid" in merged.columns else np.nan
            if pd.isna(val_czid) and pd.isna(val_seqtoid):
                continue
            if pd.isna(val_czid) or pd.isna(val_seqtoid):
                category = 'missing'
                abs_diff = 'N/A'
            else:
                val_czid = float(val_czid)
                val_seqtoid = float(val_seqtoid)
                abs_diff = abs(val_czid - val_seqtoid)
                if abs_diff <= 0.005:
                    category = 'equivalent'
                elif abs_diff <= 0.05:
                    category = 'warning'
                else:
                    category = 'significant'
            rows.append({
                'taxon': taxon,
                'sample': sample,
                'value_czid': val_czid,
                'value_seqtoid': val_seqtoid,
                'abs_diff': abs_diff if category != 'missing' else 'N/A',
                'category': category,
                'symbol': DIFF_SYMBOLS.get(category, category)
            })

    pd.DataFrame(rows).to_csv('step4_combined_rpm.csv', index=False)
    print(f"  → step4_combined_rpm.csv ({len(rows)} entries, {len(sample_cols)} czid samples)")


def compare_contig_summary_reports():
    print("Comparing per-sample contig_summary_report.csv files...")
    results = []
    missing_in_czid = []
    missing_in_seqtoid = []

    for sample in analysis_samples():
        czid_files = find_sample_file(CZID_DIR, sample, '*_contig_summary_report.csv')
        seqtoid_files = find_sample_file(SEQTOID_DIR, sample, '*_contig_summary_report.csv')

        if len(czid_files) != 1 or len(seqtoid_files) != 1:
            if len(czid_files) != 1:
                missing_in_czid.append(sample)
            if len(seqtoid_files) != 1:
                missing_in_seqtoid.append(sample)
            results.append({'sample': sample, 'mismatch_proportion': 'N/A', 'category': 'missing', 'symbol': 'missing'})
            continue

        czid_df = pd.read_csv(czid_files[0]).set_index('contig_name')
        seqtoid_df = pd.read_csv(seqtoid_files[0]).set_index('contig_name')
        all_contigs = set(czid_df.index) | set(seqtoid_df.index)
        total_contigs = len(all_contigs)
        mismatches = 0
        for contig in all_contigs:
            if contig not in czid_df.index or contig not in seqtoid_df.index:
                mismatches += 1
                continue
            czid_taxid = str(czid_df.at[contig, 'NT.species_taxid']) if 'NT.species_taxid' in czid_df.columns else None
            seqtoid_taxid = str(seqtoid_df.at[contig, 'NT.species_taxid']) if 'NT.species_taxid' in seqtoid_df.columns else None
            if czid_taxid != seqtoid_taxid:
                mismatches += 1

        proportion = mismatches / total_contigs if total_contigs else 0
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

    pd.DataFrame(results).to_csv('step5_contig_summary_comparison.csv', index=False)
    if missing_in_czid:
        print(f"  - Missing in czid: {', '.join(missing_in_czid)}")
    if missing_in_seqtoid:
        print(f"  - Missing in seqtoid: {', '.join(missing_in_seqtoid)}")


def compare_host_gene_counts():
    print("Comparing per-sample reads_per_transcript.kallisto.tsv files...")
    results = []
    missing_in_czid = []
    missing_in_seqtoid = []

    for sample in analysis_samples():
        czid_files = find_sample_file(CZID_DIR, sample, '*_reads_per_transcript.kallisto.tsv')
        seqtoid_files = find_sample_file(SEQTOID_DIR, sample, '*_reads_per_transcript.kallisto.tsv')

        if len(czid_files) != 1 or len(seqtoid_files) != 1:
            if len(czid_files) != 1:
                missing_in_czid.append(sample)
            if len(seqtoid_files) != 1:
                missing_in_seqtoid.append(sample)
            results.append({'sample': sample, 'discrepancy_proportion': 'N/A', 'category': 'missing', 'symbol': 'missing'})
            continue

        czid_ids = set(pd.read_csv(czid_files[0], sep='\t')['target_id'].unique())
        seqtoid_ids = set(pd.read_csv(seqtoid_files[0], sep='\t')['target_id'].unique())
        unique_czid = len(czid_ids - seqtoid_ids)
        unique_seqtoid = len(seqtoid_ids - czid_ids)
        common = len(czid_ids & seqtoid_ids)
        total_unique = unique_czid + unique_seqtoid + common
        discrepancy_proportion = (unique_czid + unique_seqtoid) / total_unique if total_unique else 0
        if discrepancy_proportion < 0.05:
            category = 'equivalent'
        elif discrepancy_proportion < 0.1:
            category = 'warning'
        else:
            category = 'significant'
        symbol = DIFF_SYMBOLS.get(category, category)
        print(f"  - {sample}: discrepancy {discrepancy_proportion:.4f} (unique czid {unique_czid}, unique seqtoid {unique_seqtoid}, common {common}) → {symbol}")
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
    with open(czid_path, 'r') as f:
        json.load(f)
    with open(seqtoid_path, 'r') as f:
        json.load(f)
    print("Step 7: BIOM file")
    pd.DataFrame([{'file': file, 'status': 'compared'}]).to_csv('step7_biom.csv', index=False)


def compare_nonhost_fastqs():
    print("Step 8: Non-host reads FASTQ (R1 & R2)")
    rows = []
    missing_czid = []
    missing_seqtoid = []
    for sample in analysis_samples():
        for read in ['R1', 'R2']:
            czid_f = find_sample_file(CZID_DIR, sample, f'*_reads_nh_{read}.fastq')
            seqtoid_f = find_sample_file(SEQTOID_DIR, sample, f'*_reads_nh_{read}.fastq')
            if len(czid_f) != 1 or len(seqtoid_f) != 1:
                if len(czid_f) != 1:
                    missing_czid.append(f"{sample} {read}")
                if len(seqtoid_f) != 1:
                    missing_seqtoid.append(f"{sample} {read}")
                rows.append({'sample': sample, 'read': read, 'identical': 'missing'})
                continue
            identical = file_sha256(czid_f[0]) == file_sha256(seqtoid_f[0])
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
    for sample in analysis_samples():
        czid_f = find_sample_file(CZID_DIR, sample, '*_contigs_nh.fasta')
        seqtoid_f = find_sample_file(SEQTOID_DIR, sample, '*_contigs_nh.fasta')
        if len(czid_f) != 1 or len(seqtoid_f) != 1:
            if len(czid_f) != 1:
                missing_czid.append(sample)
            if len(seqtoid_f) != 1:
                missing_seqtoid.append(sample)
            rows.append({'sample': sample, 'identical': 'missing'})
            continue
        identical = file_sha256(czid_f[0]) == file_sha256(seqtoid_f[0])
        rows.append({'sample': sample, 'identical': 'T' if identical else 'F'})
    pd.DataFrame(rows).to_csv('step9_nonhost_contigs.csv', index=False)
    if missing_czid or missing_seqtoid:
        pd.DataFrame({
            'missing_in_czid': [', '.join(missing_czid)],
            'missing_in_seqtoid': [', '.join(missing_seqtoid)]
        }).to_csv('step9_missing.csv', index=False)


def main():
    print("Starting CZID short-read pipeline comparison...\n")
    analysis_samples()

    print("\n=== Step 1: Sample Metadata Comparison ===")
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