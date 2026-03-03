#!/usr/bin/env python3
"""
CZID AMR Pipeline Comparison (updated with granular contig comparison)
"""

import os
import glob
import pandas as pd
import hashlib
from typing import List, Tuple

# ────────────────────────────────────────────────────────────────
# Configuration
# ────────────────────────────────────────────────────────────────

CZID_DIR = 'czid'
SEQTOID_DIR = 'seqtoid'

EXPECTED_SAMPLES: List[str] = [
    'ERR11417004',
    'SRR10903401',
    'SRR12876565',
    'SRR13227003',
    'SRR13227004',
    'SRR13227005',
    'SRR15049352'
]

# ────────────────────────────────────────────────────────────────
# Helpers
# ────────────────────────────────────────────────────────────────

def file_sha256(filepath: str) -> str:
    sha256 = hashlib.sha256()
    with open(filepath, 'rb') as f:
        while chunk := f.read(8192 * 1024):
            sha256.update(chunk)
    return sha256.hexdigest()


def load_contig_stats_and_seqs(path: str) -> Tuple[int, int, list]:
    """Return (num_contigs, total_bp, list_of_sequences_sorted)"""
    num = 0
    total_bp = 0
    seqs = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                num += 1
            elif line:
                seqs.append(line)
                total_bp += len(line)
    return num, total_bp, sorted(seqs)   # sorted for deterministic pairing


def compare_contigs_detailed(czid_path: str, seqtoid_path: str) -> str:
    """Returns human-readable status with granularity"""
    try:
        n1, bp1, seqs1 = load_contig_stats_and_seqs(czid_path)
        n2, bp2, seqs2 = load_contig_stats_and_seqs(seqtoid_path)

        if n1 != n2 or bp1 != bp2:
            return f"stats differ (contigs: {n1} vs {n2}, bp: {bp1} vs {bp2})"

        # same stats → check sequence content
        if file_sha256(czid_path) == file_sha256(seqtoid_path):
            return "100% identical (byte-for-byte after header strip & sort)"

        # hashes differ but stats match → compute mismatch rate (segregating sites)
        mismatches = 0
        total_aligned = 0
        for s1, s2 in zip(seqs1, seqs2):
            l = max(len(s1), len(s2))
            total_aligned += l
            if len(s1) != len(s2):
                mismatches += abs(len(s1) - len(s2))
            else:
                for a, b in zip(s1, s2):
                    if a != b:
                        mismatches += 1

        rate = mismatches / total_aligned if total_aligned > 0 else 0.0

        if rate <= 1e-8:
            return "sequences 100% identical"
        elif rate <= 0.0001:   # < 0.01%
            return f"minor differences (mismatch rate {rate:.6f})"
        elif rate <= 0.001:    # 0.01–0.1%
            return f"moderate differences (mismatch rate {rate:.6f})"
        else:
            return f"MAJOR differences (mismatch rate {rate:.6f})"

    except Exception as e:
        return f"error: {e}"


# ────────────────────────────────────────────────────────────────
# Step 1: Metadata
# ────────────────────────────────────────────────────────────────

def compare_metadata():
    file = 'sample_metadata.csv'
    czid_path = os.path.join(CZID_DIR, file)
    seqtoid_path = os.path.join(SEQTOID_DIR, file)

    if not os.path.exists(czid_path) or not os.path.exists(seqtoid_path):
        print(f"Step 1: {file} missing → skipped")
        pd.DataFrame({'file': [file], 'status': ['missing']}).to_csv('step1_metadata.csv', index=False)
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
    for sample in EXPECTED_SAMPLES:
        czid_sample = czid_df[czid_df['sample_name'] == sample]
        seqtoid_sample = seqtoid_df[seqtoid_df['sample_name'] == sample]

        czid_genes = set(czid_sample['gene_name'].dropna().unique())
        seqtoid_genes = set(seqtoid_sample['gene_name'].dropna().unique())

        missing_count = len(czid_genes - seqtoid_genes)
        total = len(czid_genes)
        proportion = missing_count / total if total > 0 else 0.0

        if proportion < 0.005:
            cat = 'equivalent'
        elif proportion < 0.05:
            cat = 'warning'
        else:
            cat = 'significant'

        results.append({
            'sample': sample,
            'total_czid_genes': total,
            'missing_in_seqtoid': missing_count,
            'missing_proportion': round(proportion, 6),
            'category': cat
        })

        print(f"  {sample}: missing proportion {proportion:.4f} ({missing_count}/{total}) → {cat}")

    pd.DataFrame(results).to_csv('step2_amr_comparison.csv', index=False)


# ────────────────────────────────────────────────────────────────
# Step 3: Non-host Contigs FASTA (with length + segregating sites)
# ────────────────────────────────────────────────────────────────

def compare_contig_fastas():
    print("\n=== Step 3: Non-host Contigs FASTA (*_contigs.fa) ===")

    rows = []
    for sample in EXPECTED_SAMPLES:
        czid_files = glob.glob(os.path.join(CZID_DIR, f"{sample}*_contigs.fa"))
        seqtoid_files = glob.glob(os.path.join(SEQTOID_DIR, f"{sample}*_contigs.fa"))

        if len(czid_files) != 1 or len(seqtoid_files) != 1:
            rows.append({'sample': sample, 'status': 'missing'})
            print(f"  {sample}: missing")
            continue

        status = compare_contigs_detailed(czid_files[0], seqtoid_files[0])
        rows.append({'sample': sample, 'status': status})
        print(f"  {sample}: {status}")

    pd.DataFrame(rows).to_csv('step3_contigs_fasta.csv', index=False)


# ────────────────────────────────────────────────────────────────
# Main
# ────────────────────────────────────────────────────────────────

def main():
    print("CZID AMR Pipeline Comparison (with granular contig check)\n")
    print(f"Samples: {', '.join(EXPECTED_SAMPLES)}\n")

    compare_metadata()
    compare_combined_amr_results()
    compare_contig_fastas()

    print("\nDone. Check the step*.csv files for full results.")


if __name__ == '__main__':
    main()