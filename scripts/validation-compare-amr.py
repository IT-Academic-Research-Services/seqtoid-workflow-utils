#!/usr/bin/env python3
"""
CZID AMR Pipeline Comparison (updated with granular contig comparison)
"""

import os
import glob
import pandas as pd
import hashlib
from typing import List, Tuple
import csv
from collections import defaultdict

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



def read_fasta_to_dict(fasta_path: str) -> dict:
    """Read FASTA → {header_first_word: sequence}"""
    d = {}
    current_id = None
    seq_lines = []
    with open(fasta_path) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if current_id:
                    d[current_id] = ''.join(seq_lines)
                current_id = line[1:].split(maxsplit=1)[0]  # first token after >
                seq_lines = []
            else:
                seq_lines.append(line)
        if current_id:
            d[current_id] = ''.join(seq_lines)
    return d


def compare_contigs_by_id_and_content(
        czid_path: str,
        seqtoid_path: str,
        sample: str,
        detail_csv: str = "step3_contigs_fasta.csv"
) -> dict:
    """
    Compare two contigs FASTA files:
      - By exact contig ID/name match
      - By exact sequence content (ignoring IDs/headers)

    Appends results to detail_csv
    Returns summary dict for console + high-level summary
    """
    czid_dict = read_fasta_to_dict(czid_path)
    seqtoid_dict = read_fasta_to_dict(seqtoid_path)

    total_czid = len(czid_dict)
    total_seqtoid = len(seqtoid_dict)

    # ─── 1. Matching by exact contig ID ───────────────────────────────
    common_ids = set(czid_dict) & set(seqtoid_dict)
    only_czid_ids = set(czid_dict) - set(seqtoid_dict)
    only_seqtoid_ids = set(seqtoid_dict) - set(czid_dict)

    snps_total = 0
    bases_total = 0
    length_mismatches = 0

    detail_rows = []

    for cid in sorted(common_ids):
        s1 = czid_dict[cid]
        s2 = seqtoid_dict[cid]
        l1, l2 = len(s1), len(s2)

        if l1 != l2:
            length_mismatches += 1
            bases_total += max(l1, l2)
            status = f"length mismatch ({l1} vs {l2})"
            snps = 0
        else:
            diffs = sum(a != b for a, b in zip(s1, s2))
            snps_total += diffs
            bases_total += l1
            status = "identical" if diffs == 0 else f"{diffs} SNPs"

        detail_rows.append([sample, "by_name", cid, status, diffs, l1])

    snp_rate_name = snps_total / bases_total if bases_total > 0 else 0.0

    by_name = {
        "total_czid": total_czid,
        "total_seqtoid": total_seqtoid,
        "common": len(common_ids),
        "only_czid": len(only_czid_ids),
        "only_seqtoid": len(only_seqtoid_ids),
        "length_mismatches": length_mismatches,
        "snps": snps_total,
        "snp_rate": round(snp_rate_name, 6),
    }

    # ─── 2. Matching by exact sequence content ────────────────────────
    czid_seq_set = set(czid_dict.values())
    seqtoid_seq_set = set(seqtoid_dict.values())
    common_sequences = czid_seq_set & seqtoid_seq_set

    by_content = {
        "common": len(common_sequences),
        "only_czid": len(czid_seq_set - seqtoid_seq_set),
        "only_seqtoid": len(seqtoid_seq_set - czid_seq_set),
    }

    # ─── Verdict logic ────────────────────────────────────────────────
    verdict = "identical"
    if by_content["only_czid"] > 0 or by_content["only_seqtoid"] > 0:
        verdict = "different contig set"
    elif by_name["snp_rate"] > 1e-5 or by_name["length_mismatches"] > 0:
        verdict = f"SNPs or length differences (rate {by_name['snp_rate']:.6f})"

    # ─── Console output ───────────────────────────────────────────────
    print(f"  {sample:20} {verdict}")
    print(f"     by name:     {by_name['common']}/{by_name['total_czid']} common "
          f"({by_name['only_czid']} only czid, {by_name['only_seqtoid']} only seqtoid) "
          f" | snp rate {by_name['snp_rate']:.6f}")
    print(f"     by content:  {by_content['common']} common sequences "
          f"({by_content['only_czid']} only czid, {by_content['only_seqtoid']} only seqtoid)")

    # ─── Append detail rows to CSV ────────────────────────────────────
    header = ["sample", "match_type", "contig_id", "status", "snps_or_gaps", "length"]
    file_exists = os.path.exists(detail_csv)

    with open(detail_csv, 'a', newline='') as f:
        writer = csv.writer(f)
        if not file_exists:
            writer.writerow(header)
        writer.writerows(detail_rows)

    # ─── Return summary for high-level CSV ────────────────────────────
    return {
        "sample": sample,
        "verdict": verdict,
        "by_name": by_name,
        "by_content": by_content
    }



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
    print("\n=== Step 3: Contig FASTA comparison ===\n")

    detail_csv = "step3_contigs_fasta.csv"
    summary_csv = "step3_contigs_summary.csv"

    # Prepare / clear detail CSV
    if os.path.exists(detail_csv):
        os.remove(detail_csv)

    summary_data = []

    for sample in EXPECTED_SAMPLES:
        cz_files = glob.glob(os.path.join(CZID_DIR, f"{sample}*_contigs.fa"))
        sq_files = glob.glob(os.path.join(SEQTOID_DIR, f"{sample}*_contigs.fa"))

        if len(cz_files) != 1 or len(sq_files) != 1:
            print(f"  {sample:20} missing or multiple files")
            summary_data.append({
                "sample": sample,
                "verdict": "missing file(s)",
                "by_name": {"total_czid": "N/A"},
                "by_content": {"common": "N/A"}
            })
            continue

        res = compare_contigs_by_id_and_content(cz_files[0], sq_files[0], sample, detail_csv)
        summary_data.append(res)

    # Write high-level summary CSV
    with open(summary_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "sample", "verdict",
            "by_name_total_czid", "by_name_common", "by_name_snp_rate",
            "by_content_common"
        ])
        writer.writeheader()
        for row in summary_data:
            writer.writerow({
                "sample": row["sample"],
                "verdict": row["verdict"],
                "by_name_total_czid": row["by_name"].get("total_czid", "N/A"),
                "by_name_common": row["by_name"].get("common", "N/A"),
                "by_name_snp_rate": row["by_name"].get("snp_rate", "N/A"),
                "by_content_common": row["by_content"].get("common", "N/A")
            })

    print(f"\nDetail   → {detail_csv}")
    print(f"Summary  → {summary_csv}")



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