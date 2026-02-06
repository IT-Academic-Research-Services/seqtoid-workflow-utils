import pandas as pd
import os
import glob
import numpy as np
import hashlib
import subprocess
from typing import List

# ────────────────────────────────────────────────────────────────
# Configuration
# ────────────────────────────────────────────────────────────────

CZID_DIR = 'czid'
SEQTOID_DIR = 'seqtoid'

EXPECTED_SAMPLES: List[str] = [
    'SRR11454628',
    'SRR34692683',
    'SRR34692681'
]

METADATA_FILE = 'sample_metadata.csv'
OVERVIEW_FILE = 'consensus_genome_overview.csv'
REPORT_FILE = 'report.tsv'
DEPTH_FILE = 'samtools_depth.txt'
CONSENSUS_PATTERN = '*_consensus.fa'  # e.g., SRR34692681_881249_MN908947.3_consensus.fa

# ────────────────────────────────────────────────────────────────
# Helpers
# ────────────────────────────────────────────────────────────────

def file_sha256(filepath):
    sha256 = hashlib.sha256()
    with open(filepath, 'rb') as f:
        while chunk := f.read(8192 * 1024):
            sha256.update(chunk)
    return sha256.hexdigest()

def numeric_diff(a: np.ndarray, b: np.ndarray, atol: float = 0.005) -> str:
    if len(a) == 0 or len(b) == 0:
        return 'empty'
    diff = np.abs(a - b)
    worst = diff.max()
    if np.isnan(worst):
        return 'NaN'
    if worst <= 0.005:
        return 'equivalent'
    elif worst <= 0.05:
        return 'warning'
    else:
        return 'significant'

def compare_numeric_dfs(df1: pd.DataFrame, df2: pd.DataFrame, id_cols: List[str]) -> dict:
    num_cols = df1.select_dtypes(include=[np.number]).columns.intersection(df2.columns)
    results = {}
    for col in num_cols:
        vals1 = df1[col].fillna(0).values.astype(float)
        vals2 = df2[col].fillna(0).values.astype(float)
        cat = numeric_diff(vals1, vals2)
        results[col] = cat
    return results

def compare_numeric_series(s1: pd.Series, s2: pd.Series) -> str:
    if len(s1) != len(s2):
        return 'length mismatch'
    vals1 = s1.fillna(0).values.astype(float)
    vals2 = s2.fillna(0).values.astype(float)
    return numeric_diff(vals1, vals2)

# ────────────────────────────────────────────────────────────────
# Step 1: Metadata
# ────────────────────────────────────────────────────────────────

def compare_metadata():
    file = METADATA_FILE
    czid_path = os.path.join(CZID_DIR, file)
    seqtoid_path = os.path.join(SEQTOID_DIR, file)

    if not os.path.exists(czid_path) or not os.path.exists(seqtoid_path):
        print(f'Error: {file} missing.')
        pd.DataFrame({'file': [file], 'identical': ['missing']}).to_csv('step1_metadata.csv', index=False)
        return

    czid_df = pd.read_csv(czid_path)
    seqtoid_df = pd.read_csv(seqtoid_path)

    sort_key = 'sample_name'
    czid_s = czid_df.sort_values(sort_key).reset_index(drop=True)
    seqtoid_s = seqtoid_df.sort_values(sort_key).reset_index(drop=True)

    identical = czid_s.equals(seqtoid_s)
    pd.DataFrame({'file': [file], 'identical': ['T' if identical else 'F']}).to_csv('step1_metadata.csv', index=False)

    if not identical:
        diffs = compare_numeric_dfs(czid_s, seqtoid_s, id_cols=[sort_key])
        if diffs:
            pd.DataFrame(list(diffs.items()), columns=['column', 'diff_category']).to_csv('step1_metadata_diffs.csv', index=False)

    # Missing/extra
    rows = []
    for name, df in [('czid', czid_df), ('seqtoid', seqtoid_df)]:
        actual = set(df['sample_name'].astype(str))
        missing = ', '.join(sorted(set(EXPECTED_SAMPLES) - actual))
        extra = ', '.join(sorted(actual - set(EXPECTED_SAMPLES)))
        if missing or extra:
            rows.append({'directory': name, 'missing': missing, 'extra': extra})
    if rows:
        pd.DataFrame(rows).to_csv('step1_missing_extra.csv', index=False)

# ────────────────────────────────────────────────────────────────
# Step 2: Overview
# ────────────────────────────────────────────────────────────────

def compare_overview():
    file = OVERVIEW_FILE
    czid_path = os.path.join(CZID_DIR, file)
    seqtoid_path = os.path.join(SEQTOID_DIR, file)

    if not os.path.exists(czid_path) or not os.path.exists(seqtoid_path):
        print(f'Error: {file} missing.')
        pd.DataFrame({'file': [file], 'identical': ['missing']}).to_csv('step2_overview.csv', index=False)
        return

    czid_df = pd.read_csv(czid_path)
    seqtoid_df = pd.read_csv(seqtoid_path)

    sort_key = 'Sample Name'
    czid_s = czid_df.sort_values(sort_key).reset_index(drop=True)
    seqtoid_s = seqtoid_df.sort_values(sort_key).reset_index(drop=True)

    identical = czid_s.equals(seqtoid_s)
    pd.DataFrame({'file': [file], 'identical': ['T' if identical else 'F']}).to_csv('step2_overview.csv', index=False)

    if not identical:
        diffs = compare_numeric_dfs(czid_s, seqtoid_s, id_cols=[sort_key])
        if diffs:
            pd.DataFrame(list(diffs.items()), columns=['column', 'diff_category']).to_csv('step2_overview_diffs.csv', index=False)

    # Missing/extra samples
    rows = []
    for name, df in [('czid', czid_df), ('seqtoid', seqtoid_df)]:
        actual = set(df['Sample Name'].astype(str))
        missing = ', '.join(sorted(set(EXPECTED_SAMPLES) - actual))
        extra = ', '.join(sorted(actual - set(EXPECTED_SAMPLES)))
        if missing or extra:
            rows.append({'directory': name, 'missing': missing, 'extra': extra})
    if rows:
        pd.DataFrame(rows).to_csv('step2_missing_extra.csv', index=False)

# ────────────────────────────────────────────────────────────────
# Step 3: report.tsv per-sample
# ────────────────────────────────────────────────────────────────

def compare_reports():
    rows = []
    missing_czid = []
    missing_seqtoid = []

    for sample in EXPECTED_SAMPLES:
        czid_sub = glob.glob(os.path.join(CZID_DIR, f'{sample}_*'))
        seqtoid_sub = glob.glob(os.path.join(SEQTOID_DIR, f'{sample}_*'))

        if len(czid_sub) != 1 or len(seqtoid_sub) != 1:
            if len(czid_sub) != 1: missing_czid.append(sample)
            if len(seqtoid_sub) != 1: missing_seqtoid.append(sample)
            rows.append({'sample': sample, 'identical': 'missing'})
            continue

        czid_path = os.path.join(czid_sub[0], REPORT_FILE)
        seqtoid_path = os.path.join(seqtoid_sub[0], REPORT_FILE)

        if not os.path.exists(czid_path) or not os.path.exists(seqtoid_path):
            rows.append({'sample': sample, 'identical': 'missing'})
            continue

        czid_df = pd.read_csv(czid_path, sep='\t', header=None, names=['Metric', 'Value'])
        seqtoid_df = pd.read_csv(seqtoid_path, sep='\t', header=None, names=['Metric', 'Value'])

        czid_s = czid_df.sort_values('Metric').reset_index(drop=True)
        seqtoid_s = seqtoid_df.sort_values('Metric').reset_index(drop=True)

        identical = czid_s.equals(seqtoid_s)
        rows.append({'sample': sample, 'identical': 'T' if identical else 'F'})

        if not identical:
            czid_df['Value_num'] = pd.to_numeric(czid_df['Value'], errors='coerce')
            seqtoid_df['Value_num'] = pd.to_numeric(seqtoid_df['Value'], errors='coerce')
            czid_s_num = czid_df.sort_values('Metric')['Value_num'].fillna(0).values
            seqtoid_s_num = seqtoid_df.sort_values('Metric')['Value_num'].fillna(0).values
            diff_cat = numeric_diff(czid_s_num, seqtoid_s_num)
            pd.DataFrame({'sample': [sample], 'overall_diff_category': [diff_cat]}).to_csv(f'step3_report_diffs_{sample}.csv', index=False)

    pd.DataFrame(rows).to_csv('step3_reports.csv', index=False)

    if missing_czid or missing_seqtoid:
        pd.DataFrame({
            'missing_in_czid': [', '.join(missing_czid)],
            'missing_in_seqtoid': [', '.join(missing_seqtoid)]
        }).to_csv('step3_missing.csv', index=False)

# ────────────────────────────────────────────────────────────────
# Step 4: FASTA
# ────────────────────────────────────────────────────────────────

def compare_fasta():
    rows = []
    missing_czid = []
    missing_seqtoid = []

    for sample in EXPECTED_SAMPLES:
        czid_f = glob.glob(os.path.join(CZID_DIR, f'{sample}{CONSENSUS_PATTERN}'))
        seqtoid_f = glob.glob(os.path.join(SEQTOID_DIR, f'{sample}{CONSENSUS_PATTERN}'))

        if len(czid_f) != 1 or len(seqtoid_f) != 1:
            if len(czid_f) != 1: missing_czid.append(sample)
            if len(seqtoid_f) != 1: missing_seqtoid.append(sample)
            rows.append({'sample': sample, 'identical': 'missing'})
            continue

        czid_path = czid_f[0]
        seqtoid_path = seqtoid_f[0]

        identical = file_sha256(czid_path) == file_sha256(seqtoid_path)
        rows.append({'sample': sample, 'identical': 'T' if identical else 'F'})

        if not identical:
            # Optional: seqkit content check
            try:
                cmd_czid = f"seqkit seq -s -i '{czid_path}' | sort | sha256sum"
                cmd_seqtoid = f"seqkit seq -s -i '{seqtoid_path}' | sort | sha256sum"
                czid_hash = subprocess.check_output(cmd_czid, shell=True, text=True).split()[0]
                seqtoid_hash = subprocess.check_output(cmd_seqtoid, shell=True, text=True).split()[0]
                content_identical = czid_hash == seqtoid_hash
                pd.DataFrame({'sample': [sample], 'content_identical': ['T' if content_identical else 'F']}).to_csv(f'step4_fasta_content_{sample}.csv', index=False)
            except:
                pass

    pd.DataFrame(rows).to_csv('step4_fasta.csv', index=False)

    if missing_czid or missing_seqtoid:
        pd.DataFrame({
            'missing_in_czid': [', '.join(missing_czid)],
            'missing_in_seqtoid': [', '.join(missing_seqtoid)]
        }).to_csv('step4_missing.csv', index=False)

# ────────────────────────────────────────────────────────────────
# Step 5: Depth
# ────────────────────────────────────────────────────────────────

def compare_depths():
    rows = []
    missing_czid = []
    missing_seqtoid = []

    for sample in EXPECTED_SAMPLES:
        czid_sub = glob.glob(os.path.join(CZID_DIR, f'{sample}_*'))
        seqtoid_sub = glob.glob(os.path.join(SEQTOID_DIR, f'{sample}_*'))

        if len(czid_sub) != 1 or len(seqtoid_sub) != 1:
            if len(czid_sub) != 1: missing_czid.append(sample)
            if len(seqtoid_sub) != 1: missing_seqtoid.append(sample)
            rows.append({'sample': sample, 'identical': 'missing'})
            continue

        czid_path = os.path.join(czid_sub[0], DEPTH_FILE)
        seqtoid_path = os.path.join(seqtoid_sub[0], DEPTH_FILE)

        if not os.path.exists(czid_path) or not os.path.exists(seqtoid_path):
            rows.append({'sample': sample, 'identical': 'missing'})
            continue

        czid_depth = pd.read_csv(czid_path, header=None, squeeze=True, dtype=float)
        seqtoid_depth = pd.read_csv(seqtoid_path, header=None, squeeze=True, dtype=float)

        if len(czid_depth) != len(seqtoid_depth):
            rows.append({'sample': sample, 'identical': 'length mismatch'})
            continue

        identical = czid_depth.equals(seqtoid_depth)
        rows.append({'sample': sample, 'identical': 'T' if identical else 'F'})

        if not identical:
            diff_cat = compare_numeric_series(czid_depth, seqtoid_depth)
            pd.DataFrame({'sample': [sample], 'overall_diff_category': [diff_cat]}).to_csv(f'step5_depth_diffs_{sample}.csv', index=False)

    pd.DataFrame(rows).to_csv('step5_depths.csv', index=False)

    if missing_czid or missing_seqtoid:
        pd.DataFrame({
            'missing_in_czid': [', '.join(missing_czid)],
            'missing_in_seqtoid': [', '.join(missing_seqtoid)]
        }).to_csv('step5_missing.csv', index=False)

# ────────────────────────────────────────────────────────────────
# Main
# ────────────────────────────────────────────────────────────────

def main():
    print('Starting CZID consensus genome pipeline comparison...\n')

    print('=== Step 1: Sample Metadata Comparison ===')
    compare_metadata()

    print('\n=== Step 2: Consensus Genome Overview Comparison ===')
    compare_overview()

    print('\n=== Step 3: Per-sample Report.tsv Comparison ===')
    compare_reports()

    print('\n=== Step 4: Consensus FASTA Comparison ===')
    compare_fasta()

    print('\n=== Step 5: Per-sample Samtools Depth Comparison ===')
    compare_depths()

    print('\nComparison complete. Check CSV files.')

if __name__ == '__main__':
    main()