import os
import argparse
import json
import statistics
from pysam import VariantFile
from Bio import SeqIO
from typing import List


def find_rust_variants_file(directory):
    for file in os.listdir(directory):
        if file.endswith('_variants.bcf'):
            return os.path.join(directory, file)
    raise FileNotFoundError("No variants.bcf file found in the Rust directory")


def find_wdl_variants_file(directory):
    path = os.path.join(directory, 'out/call_variants_out_variants_ch/variants.vcf.gz')
    if os.path.exists(path):
        return path
    raise FileNotFoundError("No variants.vcf.gz file found in the WDL directory")


def load_variants(file_path):
    variants = {}
    with VariantFile(file_path) as vcf:
        for rec in vcf:
            contig = rec.contig
            pos = rec.pos
            ref = rec.ref
            alts = tuple(sorted(rec.alts))  # Sort for consistent comparison
            key = (contig, pos)
            if key in variants:
                raise ValueError(f"Duplicate position {key} found in {file_path}")
            variants[key] = (ref, alts)
    return variants


def compare_variants(rust_vars, wdl_vars):
    rust_keys = set(rust_vars.keys())
    wdl_keys = set(wdl_vars.keys())

    common = rust_keys & wdl_keys
    rust_only = rust_keys - wdl_keys
    wdl_only = wdl_keys - rust_keys

    mismatches = []
    for key in common:
        if rust_vars[key] != wdl_vars[key]:
            mismatches.append(key)

    print("=== Variants Comparison ===")
    print(f"Common positions: {len(common)}")
    print(f"Positions unique to Rust: {len(rust_only)}")
    if rust_only:
        print("Rust-only positions:")
        for key in sorted(rust_only):
            print(f"  {key}: REF={rust_vars[key][0]}, ALT={','.join(rust_vars[key][1])}")

    print(f"Positions unique to WDL: {len(wdl_only)}")
    if wdl_only:
        print("WDL-only positions:")
        for key in sorted(wdl_only):
            print(f"  {key}: REF={wdl_vars[key][0]}, ALT={','.join(wdl_vars[key][1])}")

    print(f"Mismatches in common positions: {len(mismatches)}")
    if mismatches:
        print("Mismatches:")
        for key in sorted(mismatches):
            print(f"  {key}:")
            print(f"    Rust: REF={rust_vars[key][0]}, ALT={','.join(rust_vars[key][1])}")
            print(f"    WDL: REF={wdl_vars[key][0]}, ALT={','.join(wdl_vars[key][1])}")
    print()


def find_rust_consensus_file(directory):
    for file in os.listdir(directory):
        if file.endswith('_consensus.fa'):
            return os.path.join(directory, file)
    raise FileNotFoundError("No consensus.fa file found in the Rust directory")


def find_wdl_consensus_file(directory):
    path = os.path.join(directory, 'out/make_consensus_out_consensus_fa/consensus.fa')
    if os.path.exists(path):
        return path
    raise FileNotFoundError("No consensus.fa file found in the WDL directory")


def load_fasta(file_path):
    records = list(SeqIO.parse(file_path, "fasta"))
    if len(records) == 0:
        raise ValueError(f"No sequences found in {file_path}")
    if len(records) > 1:
        raise ValueError(f"Multiple sequences found in {file_path}; expected exactly one")
    return str(records[0].seq).upper()


def compare_consensus(rust_dir, wdl_dir):
    rust_file = find_rust_consensus_file(rust_dir)
    wdl_file = find_wdl_consensus_file(wdl_dir)

    print("=== Consensus Comparison ===")
    rust_seq = load_fasta(rust_file)
    wdl_seq = load_fasta(wdl_file)

    print(f"Rust sequence length: {len(rust_seq)}")
    print(f"WDL sequence length: {len(wdl_seq)}")

    if rust_seq == wdl_seq:
        print("Sequences are identical")
    else:
        print("Sequences differ")
        for i, (r, w) in enumerate(zip(rust_seq, wdl_seq)):
            if r != w:
                print(f"  Position {i + 1}: Rust={r}, WDL={w}")
                break  # Show only first mismatch for brevity
        if len(rust_seq) != len(wdl_seq):
            print(f"  Length mismatch: Rust={len(rust_seq)}, WDL={len(wdl_seq)}")
    print()


def find_rust_stats_file(directory):
    for file in os.listdir(directory):
        if file.endswith('_stats.json'):
            return os.path.join(directory, file)
    raise FileNotFoundError("No stats.json file found in the Rust directory")


def find_wdl_stats_file(directory):
    path = os.path.join(directory, 'out/compute_stats_out_output_stats/stats.json')
    if os.path.exists(path):
        return path
    raise FileNotFoundError("No stats.json file found in the WDL directory")


def find_wdl_depth_file(directory):
    path = os.path.join(directory, 'out/compute_stats_out_sam_depths/samtools_depth.txt')
    if os.path.exists(path):
        return path
    raise FileNotFoundError("No samtools_depth.txt file found in the WDL directory")


def load_stats(file_path, is_wdl=False):
    with open(file_path, 'r') as f:
        stats = json.load(f)
    if is_wdl:
        # Map WDL dotted keys to Rust-style keys
        mapping = {
            'depth_q.25': 'depth_q25',
            'depth_q.5': 'depth_q50',
            'depth_q.75': 'depth_q75',
        }
        for old_key, new_key in mapping.items():
            if old_key in stats:
                stats[new_key] = stats.pop(old_key)
    return stats


def load_depth_file(file_path) -> List[int]:
    depths = []
    with open(file_path, 'r') as f:
        for line in f:
            stripped_line = line.strip()
            if stripped_line:
                try:
                    depth = int(stripped_line)
                    depths.append(depth)
                except ValueError:
                    raise ValueError(f"Invalid depth value in {file_path}: {stripped_line}")
    return depths


def compute_depth_stats(depths: List[int]):
    if not depths:
        return {
            'depth_avg': 0.0,
            'depth_q25': 0.0,
            'depth_q50': 0.0,
            'depth_q75': 0.0,
            'depth_frac_above_10x': 0.0,
            'depth_frac_above_25x': 0.0,
            'depth_frac_above_50x': 0.0,
            'depth_frac_above_100x': 0.0
        }
    depth_avg = statistics.mean(depths)
    depth_q25, depth_q50, depth_q75 = statistics.quantiles(depths, n=4)
    total_positions = len(depths)
    depth_frac_above_10x = sum(1 for d in depths if d >= 10) / total_positions
    depth_frac_above_25x = sum(1 for d in depths if d >= 25) / total_positions
    depth_frac_above_50x = sum(1 for d in depths if d >= 50) / total_positions
    depth_frac_above_100x = sum(1 for d in depths if d >= 100) / total_positions
    return {
        'depth_avg': depth_avg,
        'depth_q25': depth_q25,
        'depth_q50': depth_q50,
        'depth_q75': depth_q75,
        'depth_frac_above_10x': depth_frac_above_10x,
        'depth_frac_above_25x': depth_frac_above_25x,
        'depth_frac_above_50x': depth_frac_above_50x,
        'depth_frac_above_100x': depth_frac_above_100x
    }


def compare_stats(rust_dir, wdl_dir):
    rust_file = find_rust_stats_file(rust_dir)
    wdl_file = find_wdl_stats_file(wdl_dir)
    wdl_depth_file = find_wdl_depth_file(wdl_dir)

    print("=== Stats Comparison ===")
    rust_stats = load_stats(rust_file)
    wdl_stats = load_stats(wdl_file, is_wdl=True)  # Apply mapping for WDL
    wdl_depths = load_depth_file(wdl_depth_file)

    # Compute depth statistics from WDL depth file
    wdl_depth_stats = compute_depth_stats(wdl_depths)

    # Fields to compare (including depth stats from Rust stats.json)
    fields = [
        'depth_avg', 'depth_q25', 'depth_q50', 'depth_q75',
        'depth_frac_above_10x', 'depth_frac_above_25x', 'depth_frac_above_50x', 'depth_frac_above_100x',
        'total_reads', 'mapped_reads', 'mapped_paired', 'paired_inward', 'paired_outward',
        'paired_other_orientation', 'ercc_mapped_reads', 'ercc_mapped_paired',
        'ref_snps', 'ref_mnps', 'ref_indels', 'n_actg', 'n_missing', 'n_gap', 'n_ambiguous',
        'coverage_breadth', 'max_aligned_length', 'total_length', 'coverage_bin_size'
    ]

    print("Comparing key statistics:")
    mismatches = []
    for field in fields:
        rust_value = rust_stats.get(field)
        wdl_value = wdl_stats.get(field, wdl_depth_stats.get(field))

        # Handle None values (e.g., optional fields like paired_inward)
        if rust_value is None and wdl_value is None:
            continue
        elif rust_value is None or wdl_value is None:
            mismatches.append((field, rust_value, wdl_value))
            continue

        # Compare numeric values with tolerance for floating-point
        if isinstance(rust_value, (int, float)) and isinstance(wdl_value, (int, float)):
            if abs(rust_value - wdl_value) > 1e-6:  # Tolerance for floating-point
                mismatches.append((field, rust_value, wdl_value))
        elif rust_value != wdl_value:
            mismatches.append((field, rust_value, wdl_value))

    # Compare allele_counts separately (dictionary comparison)
    rust_allele_counts = rust_stats.get('allele_counts', {})
    wdl_allele_counts = wdl_stats.get('allele_counts', {})
    allele_mismatches = []
    allele_keys = set(rust_allele_counts.keys()) | set(wdl_allele_counts.keys())
    for key in allele_keys:
        rust_count = rust_allele_counts.get(key, 0)
        wdl_count = wdl_allele_counts.get(key, 0)
        if rust_count != wdl_count:
            allele_mismatches.append((key, rust_count, wdl_count))

    # Report results
    print(f"Fields compared: {len(fields)}")
    print(f"Mismatches in fields: {len(mismatches)}")
    if mismatches:
        print("Field mismatches:")
        for field, rust_val, wdl_val in mismatches:
            print(f"  {field}: Rust={rust_val}, WDL={wdl_val}")

    print(f"Allele count mismatches: {len(allele_mismatches)}")
    if allele_mismatches:
        print("Allele count mismatches:")
        for key, rust_count, wdl_count in sorted(allele_mismatches):
            print(f"  {key}: Rust={rust_count}, WDL={wdl_count}")

    print("Note: 'coverage' array comparison skipped (requires specific handling)")
    print("Note: WDL quantile keys were mapped from 'depth_q.XX' to 'depth_qXX' for comparison")
    print()


def validate_pipeline(rust_dir, wdl_dir):
    # Run all comparisons
    rust_vars = load_variants(find_rust_variants_file(rust_dir))
    wdl_vars = load_variants(find_wdl_variants_file(wdl_dir))
    compare_variants(rust_vars, wdl_vars)

    compare_consensus(rust_dir, wdl_dir)
    compare_stats(rust_dir, wdl_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Validate Rust pipeline results against WDL pipeline")
    parser.add_argument('rust_dir', help='Top-level result directory for Rust pipeline')
    parser.add_argument('wdl_dir', help='Top-level result directory for WDL pipeline')

    args = parser.parse_args()
    validate_pipeline(args.rust_dir, args.wdl_dir)