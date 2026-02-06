import pandas as pd
import os
import glob
import hashlib
import subprocess
from typing import List

# ────────────────────────────────────────────────────────────────
# Configuration
# ────────────────────────────────────────────────────────────────

CZID_DIR = 'czid'
SEQTOID_DIR = 'seqtoid'

# Explicit list — no metadata loading
EXPECTED_SAMPLES: List[str] = [
    'ERR11417004',
    'SRR10903401',
    'SRR12876565',
    'SRR13227003',
    'SRR13227004',
    'SRR13227005',
    'SRR15049352'
    # add any others you need, e.g.:
    # 'SRR23038836_75M_1',
    # 'SRR11454628',
    # ...
]

DIFF_SYMBOLS = {
    'equivalent':   '✅ <0.005',
    'warning':      '⚠️ 0.005–0.05',
    'significant':  '❌ >0.05',
    'identical':    'T',
    'differ':       'F',
    'missing':      'missing',
    'content_identical': 'sequences identical',
    'content_diff': 'sequences differ',
}

# ────────────────────────────────────────────────────────────────
# Helpers
# ────────────────────────────────────────────────────────────────

def file_sha256(filepath: str) -> str:
    sha256 = hashlib.sha256()
    with open(filepath, 'rb') as f:
        while chunk := f.read(8192 * 1024):
            sha256.update(chunk)
    return sha256.hexdigest()


# ────────────────────────────────────────────────────────────────
# Step 1: Metadata (optional but kept for consistency)
# ────────────────────────────────────────────────────────────────

def compare_metadata():
    file = 'sample_metadata.csv'
    czid_path = os.path.join(CZID_DIR, file)
    seqtoid_path = os.path.join(SEQTOID_DIR, file)

    if not os.path.exists(czid_path) or not os.path.exists(seqtoid_path):
        print(f"Step 1: {file} missing in one or both dirs → skipped")
        pd.DataFrame({'file': [file], 'identical': ['missing']}).to_csv('step1_metadata.csv', index=False)
        return

    czid_df = pd.read_csv(czid_path)
    seqtoid_df = pd.read_csv(seqtoid_path)

    czid_s = czid_df.sort_values('sample_name').reset_index(drop=True)
    seqtoid_s = seqtoid_df.sort_values('sample_name').reset_index(drop=True)

    identical = czid_s.equals(seqtoid_s)
    pd.DataFrame({'file': [file], 'identical': ['T' if identical else 'F']}).to_csv('step1_metadata.csv', index=False)

    print(f"Step 1: metadata identical = {'T' if identical else 'F'}")


# ────────────────────────────────────────────────────────────────
# Step 2: Combined AMR Results (gene presence)
# ────────────────────────────────────────────────────────────────

def compare_combined_amr_results():
    file_name = 'combined_amr_results.csv'
    czid_path = os.path.join(CZID_DIR, file_name)
    seqtoid_path = os.path.join(SEQTOID_DIR, file_name)

    print("\n=== Step 2: Combined AMR Results (gene presence) ===")

    if not os.path.exists(czid_path) or not os.path.exists(seqtoid_path):
        print(f"  → {file_name} missing")
        pd.DataFrame({'file': [file_name], 'status': ['missing']}).to_csv('step2_amr_comparison.csv', index=False)
        return

    czid_df = pd.read_csv(czid_path)
    seqtoid_df = pd.read_csv(seqtoid_path)

    results = []
    missing_czid = []
    missing_seqtoid = []

    for sample in EXPECTED_SAMPLES:
        czid_sample = czid_df[czid_df['sample_name'] == sample]
        seqtoid_sample = seqtoid_df[seqtoid_df['sample_name'] == sample]

        if czid_sample.empty:
            missing_czid.append(sample)
        if seqtoid_sample.empty:
            missing_seqtoid.append(sample)

        if czid_sample.empty or seqtoid_sample.empty:
            results.append({'sample': sample, 'missing_proportion': 'N/A', 'category': 'missing', 'symbol': 'missing'})
            continue

        czid_genes = set(czid_sample['gene_name'].dropna().unique())
        seqtoid_genes = set(seqtoid_sample['gene_name'].dropna().unique())

        missing_in_seqtoid_count = len(czid_genes - seqtoid_genes)
        total_czid = len(czid_genes)
        proportion = missing_in_seqtoid_count / total_czid if total_czid > 0 else 0.0

        if proportion < 0.005:
            cat = 'equivalent'
        elif proportion < 0.05:
            cat = 'warning'
        else:
            cat = 'significant'

        symbol = DIFF_SYMBOLS.get(cat, cat)

        results.append({
            'sample': sample,
            'total_czid_genes': total_czid,
            'missing_in_seqtoid': missing_in_seqtoid_count,
            'missing_proportion': round(proportion, 6),
            'category': cat,
            'symbol': symbol
        })

        print(f"  {sample}: missing proportion {proportion:.4f} ({missing_in_seqtoid_count}/{total_czid}) → {symbol}")

    pd.DataFrame(results).to_csv('step2_amr_comparison.csv', index=False)

    if missing_czid or missing_seqtoid:
        pd.DataFrame({
            'missing_in_czid': [', '.join(missing_czid)],
            'missing_in_seqtoid': [', '.join(missing_seqtoid)]
        }).to_csv('step2_missing.csv', index=False)


# ────────────────────────────────────────────────────────────────
# Step 3: Non-host Contigs FASTA
# ────────────────────────────────────────────────────────────────

def compare_contig_fastas():
    print("\n=== Step 3: Non-host Contigs FASTA (*_contigs.fa) ===")

    rows = []
    missing_czid = []
    missing_seqtoid = []

    for sample in EXPECTED_SAMPLES:
        czid_files = glob.glob(os.path.join(CZID_DIR, f"{sample}*_contigs.fa"))
        seqtoid_files = glob.glob(os.path.join(SEQTOID_DIR, f"{sample}*_contigs.fa"))

        if len(czid_files) != 1 or len(seqtoid_files) != 1:
            if len(czid_files) != 1: missing_czid.append(sample)
            if len(seqtoid_files) != 1: missing_seqtoid.append(sample)
            rows.append({'sample': sample, 'identical': 'missing', 'note': ''})
            print(f"  {sample}: missing")
            continue

        czid_path = czid_files[0]
        seqtoid_path = seqtoid_files[0]

        file_identical = file_sha256(czid_path) == file_sha256(seqtoid_path)

        if file_identical:
            status = 'identical'
            note = ''
        else:
            try:
                cmd = (
                    f"seqkit seq -s '{czid_path}' | sort | sha256sum | cut -d' ' -f1 && "
                    f"seqkit seq -s '{seqtoid_path}' | sort | sha256sum | cut -d' ' -f1"
                )
                result = subprocess.check_output(cmd, shell=True, text=True).strip().split('\n')
                if len(result) == 2 and result[0] == result[1]:
                    status = 'content_identical'
                    note = '(sequences match, headers/order differ)'
                else:
                    status = 'content_diff'
                    note = '(sequences differ)'
            except:
                status = 'differ'
                note = '(content check failed)'

        rows.append({
            'sample': sample,
            'identical': DIFF_SYMBOLS.get(status, status),
            'note': note
        })

        print(f"  {sample}: {DIFF_SYMBOLS.get(status, status)} {note}")

    pd.DataFrame(rows).to_csv('step3_contigs_fasta.csv', index=False)

    if missing_czid or missing_seqtoid:
        pd.DataFrame({
            'missing_in_czid': [', '.join(missing_czid)],
            'missing_in_seqtoid': [', '.join(missing_seqtoid)]
        }).to_csv('step3_missing.csv', index=False)


# ────────────────────────────────────────────────────────────────
# Main
# ────────────────────────────────────────────────────────────────

def main():
    print("Starting CZID AMR pipeline comparison...\n")
    print(f"Samples being checked: {len(EXPECTED_SAMPLES)}")
    print(", ".join(EXPECTED_SAMPLES))
    print()

    compare_metadata()
    compare_combined_amr_results()
    compare_contig_fastas()

    print("\nDone. Check the step*.csv files.")


if __name__ == '__main__':
    main()