import os
import argparse
import json
import statistics
from pysam import VariantFile
from Bio import SeqIO
from typing import List, Dict


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
            alts = tuple(sorted(rec.alts))
            key = (contig, pos)
            if key not in variants:
                variants[key] = []
            variants[key].append((ref, alts))
    return variants


def compare_variants(rust_dir, wdl_dir):
    rust_vars = load_variants(find_rust_variants_file(rust_dir))
    wdl_vars = load_variants(find_wdl_variants_file(wdl_dir))

    rust_keys = set(rust_vars.keys())
    wdl_keys = set(wdl_vars.keys())

    common = rust_keys & wdl_keys
    rust_only = rust_keys - wdl_keys
    wdl_only = wdl_keys - rust_keys

    mismatches = []
    for key in sorted(common):
        rust_variants = set(rust_vars[key])
        wdl_variants = set(wdl_vars[key])
        if rust_variants != wdl_variants:
            mismatches.append((key, list(rust_variants), list(wdl_variants)))

    return {
        'common_positions': len(common),
        'rust_only_positions': len(rust_only),
        'wdl_only_positions': len(wdl_only),
        'mismatches': len(mismatches),
        'mismatch_details': mismatches,
        'rust_only_details': [(key, list(rust_vars[key])) for key in sorted(rust_only)],
        'wdl_only_details': [(key, list(wdl_vars[key])) for key in sorted(wdl_only)]
    }


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

    rust_seq = load_fasta(rust_file)
    wdl_seq = load_fasta(wdl_file)

    rust_length = len(rust_seq)
    wdl_length = len(wdl_seq)
    identical = rust_seq == wdl_seq
    mismatch_details = []
    if not identical:
        for i, (r, w) in enumerate(zip(rust_seq, wdl_seq)):
            if r != w:
                mismatch_details.append((i + 1, r, w))
                break

    return {
        'rust_length': rust_length,
        'wdl_length': wdl_length,
        'identical': identical,
        'mismatch_position': mismatch_details[0][0] if mismatch_details else None,
        'mismatch_details': mismatch_details
    }


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

    rust_stats = load_stats(rust_file)
    wdl_stats = load_stats(wdl_file, is_wdl=True)
    wdl_depths = load_depth_file(wdl_depth_file)

    # Adjust Rust stats for concatenated paired reads (double alignment metrics only)
    alignment_fields = ['mapped_reads', 'ercc_mapped_reads', 'ercc_mapped_paired']
    for field in alignment_fields:
        if rust_stats.get(field) is not None:
            rust_stats[field] *= 2

    # Compute depth statistics from WDL depth file
    wdl_depth_stats = compute_depth_stats(wdl_depths)

    # Fields to compare (excluding paired-end metrics unavailable in Rust)
    fields = [
        'depth_avg', 'depth_q25', 'depth_q50', 'depth_q75',
        'depth_frac_above_10x', 'depth_frac_above_25x', 'depth_frac_above_50x', 'depth_frac_above_100x',
        'total_reads', 'mapped_reads', 'ercc_mapped_reads', 'ercc_mapped_paired',
        'ref_snps', 'ref_mnps', 'ref_indels', 'n_actg', 'n_missing', 'n_gap', 'n_ambiguous',
        'coverage_breadth', 'max_aligned_length', 'total_length', 'coverage_bin_size'
    ]

    within_1_percent = []
    within_5_percent = []
    outside_5_percent = []
    mismatches = []
    for field in sorted(fields):
        rust_value = rust_stats.get(field)
        wdl_value = wdl_stats.get(field, wdl_depth_stats.get(field))

        # Handle None values
        if rust_value is None and wdl_value is None:
            continue
        elif rust_value is None or wdl_value is None:
            outside_5_percent.append((field, rust_value, wdl_value, None))
            mismatches.append((field, rust_value, wdl_value, None))
            continue

        # Compare numeric values
        if isinstance(rust_value, (int, float)) and isinstance(wdl_value, (int, float)):
            # For integer fields (counts), allow small absolute difference
            if field in ['total_reads', 'mapped_reads', 'ercc_mapped_reads', 'ercc_mapped_paired',
                         'ref_snps', 'ref_mnps', 'ref_indels', 'n_actg', 'n_missing', 'n_gap', 'n_ambiguous',
                         'max_aligned_length', 'total_length']:
                absolute_diff = abs(rust_value - wdl_value)
                if absolute_diff > 50:
                    outside_5_percent.append((field, rust_value, wdl_value, None))
                    mismatches.append((field, rust_value, wdl_value, None))
                else:
                    within_1_percent.append((field, rust_value, wdl_value, absolute_diff))
            else:
                # For floating-point fields, check relative difference
                if wdl_value != 0:
                    relative_diff = abs(rust_value - wdl_value) / abs(wdl_value) * 100
                    if relative_diff <= 1.0:
                        within_1_percent.append((field, rust_value, wdl_value, relative_diff))
                    elif relative_diff <= 5.0:
                        within_5_percent.append((field, rust_value, wdl_value, relative_diff))
                    else:
                        outside_5_percent.append((field, rust_value, wdl_value, relative_diff))
                        mismatches.append((field, rust_value, wdl_value, relative_diff))
                else:
                    if abs(rust_value) > 1e-6:
                        outside_5_percent.append((field, rust_value, wdl_value, float('inf')))
                        mismatches.append((field, rust_value, wdl_value, float('inf')))
                    else:
                        within_1_percent.append((field, rust_value, wdl_value, 0.0))

        else:
            if rust_value != wdl_value:
                outside_5_percent.append((field, rust_value, wdl_value, None))
                mismatches.append((field, rust_value, wdl_value, None))

    # Compare allele_counts separately
    rust_allele_counts = rust_stats.get('allele_counts', {})
    wdl_allele_counts = wdl_stats.get('allele_counts', {})
    allele_mismatches = []
    allele_keys = sorted(set(rust_allele_counts.keys()) | set(wdl_allele_counts.keys()))
    for key in allele_keys:
        rust_count = rust_allele_counts.get(key, 0)
        wdl_count = wdl_allele_counts.get(key, 0)
        absolute_diff = abs(rust_count - wdl_count)
        if absolute_diff > 0:  # Report all mismatches, as per previous script
            allele_mismatches.append((key, rust_count, wdl_count, absolute_diff))

    return {
        'fields_compared': len(fields),
        'within_1_percent': len(within_1_percent),
        'within_1_percent_details': within_1_percent,
        'within_5_percent': len(within_5_percent),
        'within_5_percent_details': within_5_percent,
        'outside_5_percent': len(outside_5_percent),
        'outside_5_percent_details': outside_5_percent,
        'allele_mismatches': len(allele_mismatches),
        'allele_mismatch_details': allele_mismatches
    }


def find_rust_depth_plot_file(directory):
    for file in os.listdir(directory):
        if file.endswith('_depth.png'):
            return os.path.join(directory, file)
    raise FileNotFoundError("No depth.png file found in the Rust directory")


def find_wdl_depth_plot_file(directory):
    path = os.path.join(directory, 'out/compute_stats_out_depths_fig/depths.png')
    if os.path.exists(path):
        return path
    raise FileNotFoundError("No depths.png file found in the WDL directory")


def compare_depth_plots(rust_dir, wdl_dir):
    try:
        rust_file = find_rust_depth_plot_file(rust_dir)
        wdl_file = find_wdl_depth_plot_file(wdl_dir)

        rust_size = os.path.getsize(rust_file)
        wdl_size = os.path.getsize(wdl_file)

        identical = rust_size == wdl_size
        return {
            'identical': identical,
            'rust_size': rust_size,
            'wdl_size': wdl_size
        }
    except FileNotFoundError:
        return {
            'identical': False,
            'rust_size': None,
            'wdl_size': None
        }


def find_rust_quast_report_file(directory):
    path = os.path.join(directory, 'quast/report.txt')
    if os.path.exists(path):
        return path
    raise FileNotFoundError("No quast/report.txt file found in the Rust directory")


def find_wdl_quast_report_file(directory):
    path = os.path.join(directory, 'out/quast_out_quast_txt/report.txt')
    if os.path.exists(path):
        return path
    raise FileNotFoundError("No quast_out_quast_txt/report.txt file found in the WDL directory")


def load_quast_report(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
    report = {}
    for line in lines:
        if '=' in line:
            key, value = line.split('=', 1)
            report[key.strip()] = value.strip()
    return report


def compare_quast_reports(rust_dir, wdl_dir):
    try:
        rust_file = find_rust_quast_report_file(rust_dir)
        wdl_file = find_wdl_quast_report_file(wdl_dir)

        rust_report = load_quast_report(rust_file)
        wdl_report = load_quast_report(wdl_file)

        common_keys = set(rust_report.keys()) & set(wdl_report.keys())
        mismatches = []
        for key in sorted(common_keys):
            if rust_report[key] != wdl_report[key]:
                mismatches.append((key, rust_report[key], wdl_report[key]))

        return {
            'common_keys': len(common_keys),
            'mismatches': len(mismatches),
            'mismatch_details': mismatches
        }
    except FileNotFoundError:
        return {
            'common_keys': 0,
            'mismatches': 0,
            'mismatch_details': []
        }


def find_rust_kraken_report_file(directory):
    for file in os.listdir(directory):
        if file.endswith('_kraken2_report.txt'):
            return os.path.join(directory, file)
    raise FileNotFoundError("No kraken2_report.txt file found in the Rust directory")


def find_wdl_kraken_report_file(directory):
    path = os.path.join(directory, 'out/filter_reads_out_kraken2_report/kraken2_report.txt')
    if os.path.exists(path):
        return path
    raise FileNotFoundError("No kraken2_report.txt file found in the WDL directory")


def load_kraken_report(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
    report = {}
    for line in lines:
        fields = line.strip().split()
        if len(fields) >= 6:
            try:
                report[fields[5]] = int(fields[1])
            except ValueError:
                continue
    return report


def compare_kraken_reports(rust_dir, wdl_dir):
    try:
        rust_file = find_rust_kraken_report_file(rust_dir)
        wdl_file = find_wdl_kraken_report_file(wdl_dir)

        rust_report = load_kraken_report(rust_file)
        wdl_report = load_kraken_report(wdl_file)

        common_keys = set(rust_report.keys()) & set(wdl_report.keys())
        mismatches = []
        for key in sorted(common_keys):
            if rust_report[key] != wdl_report[key]:
                mismatches.append((key, rust_report[key], wdl_report[key]))

        return {
            'common_keys': len(common_keys),
            'mismatches': len(mismatches),
            'mismatch_details': mismatches
        }
    except FileNotFoundError:
        return {
            'common_keys': 0,
            'mismatches': 0,
            'mismatch_details': []
        }


def find_rust_ercc_stats_file(directory):
    for file in os.listdir(directory):
        if file.endswith('_ercc_stats.txt'):
            return os.path.join(directory, file)
    raise FileNotFoundError("No ercc_stats.txt file found in the Rust directory")


def find_wdl_ercc_stats_file(directory):
    path = os.path.join(directory, 'out/quantify_erccs_out_ercc_out/ercc_stats.txt')
    if os.path.exists(path):
        return path
    raise FileNotFoundError("No ercc_stats.txt file found in the WDL directory")


def load_ercc_stats(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
    stats = {}
    for line in lines:
        fields = line.strip().split()
        if len(fields) >= 2:
            try:
                stats[fields[0]] = int(fields[1])
            except ValueError:
                continue
    return stats


def compare_ercc_stats(rust_dir, wdl_dir):
    try:
        rust_file = find_rust_ercc_stats_file(rust_dir)
        wdl_file = find_wdl_ercc_stats_file(wdl_dir)

        rust_stats = load_ercc_stats(rust_file)
        wdl_stats = load_ercc_stats(wdl_file)

        common_keys = set(rust_stats.keys()) & set(wdl_stats.keys())
        mismatches = []
        for key in sorted(common_keys):
            if rust_stats[key] != wdl_stats[key]:
                mismatches.append((key, rust_stats[key], wdl_stats[key]))

        return {
            'common_keys': len(common_keys),
            'mismatches': len(mismatches),
            'mismatch_details': mismatches
        }
    except FileNotFoundError:
        return {
            'common_keys': 0,
            'mismatches': 0,
            'mismatch_details': []
        }


def validate_pipeline(input_file):
    with open(input_file, 'r') as f:
        lines = f.readlines()

    all_results = []
    for line in lines:
        parts = line.strip().split()
        if len(parts) != 2:
            print(f"Skipping invalid line: {line.strip()}")
            continue
        rust_dir, wdl_dir = parts
        sample_name = os.path.basename(rust_dir.rstrip('/'))

        try:
            result = {}
            result['sample_name'] = sample_name

            try:
                result['variants'] = compare_variants(rust_dir, wdl_dir)
            except Exception as e:
                result['variants'] = {'error': str(e)}

            try:
                result['consensus'] = compare_consensus(rust_dir, wdl_dir)
            except Exception as e:
                result['consensus'] = {'error': str(e)}

            try:
                result['stats'] = compare_stats(rust_dir, wdl_dir)
            except Exception as e:
                result['stats'] = {'error': str(e)}

            result['depth_plots'] = compare_depth_plots(rust_dir, wdl_dir)
            result['quast_reports'] = compare_quast_reports(rust_dir, wdl_dir)
            result['kraken_reports'] = compare_kraken_reports(rust_dir, wdl_dir)
            result['ercc_stats'] = compare_ercc_stats(rust_dir, wdl_dir)

            with open(f"{sample_name}-validation.txt", 'w') as f:
                f.write("=== Variants Comparison ===\n")
                if 'error' in result['variants']:
                    f.write(f"Error: {result['variants']['error']}\n")
                else:
                    f.write(f"Common positions: {result['variants']['common_positions']}\n")
                    f.write(f"Positions unique to Rust: {result['variants']['rust_only_positions']}\n")
                    if result['variants']['rust_only_details']:
                        f.write("Rust-only positions:\n")
                        for key, variants in result['variants']['rust_only_details']:
                            for ref, alts in variants:
                                f.write(f"  {key}: REF={ref}, ALT={','.join(alts)}\n")
                    f.write(f"Positions unique to WDL: {result['variants']['wdl_only_positions']}\n")
                    if result['variants']['wdl_only_details']:
                        f.write("WDL-only positions:\n")
                        for key, variants in result['variants']['wdl_only_details']:
                            for ref, alts in variants:
                                f.write(f"  {key}: REF={ref}, ALT={','.join(alts)}\n")
                    f.write(f"Mismatches in common positions: {result['variants']['mismatches']}\n")
                    if result['variants']['mismatch_details']:
                        f.write("Mismatches:\n")
                        for key, rust_vars, wdl_vars in result['variants']['mismatch_details']:
                            f.write(f"  {key}:\n")
                            f.write(
                                f"    Rust: {', '.join(f'REF={ref}, ALT={','.join(alts)}' for ref, alts in rust_vars)}\n")
                            f.write(
                                f"    WDL: {', '.join(f'REF={ref}, ALT={','.join(alts)}' for ref, alts in wdl_vars)}\n")
                f.write("\n")

                f.write("=== Consensus Comparison ===\n")
                if 'error' in result['consensus']:
                    f.write(f"Error: {result['consensus']['error']}\n")
                else:
                    f.write(f"Rust sequence length: {result['consensus']['rust_length']}\n")
                    f.write(f"WDL sequence length: {result['consensus']['wdl_length']}\n")
                    if result['consensus']['identical']:
                        f.write("Sequences are identical\n")
                    else:
                        f.write("Sequences differ\n")
                        if result['consensus']['mismatch_details']:
                            pos, rust_base, wdl_base = result['consensus']['mismatch_details'][0]
                            f.write(f"  Position {pos}: Rust={rust_base}, WDL={wdl_base}\n")
                        if result['consensus']['rust_length'] != result['consensus']['wdl_length']:
                            f.write(
                                f"  Length mismatch: Rust={result['consensus']['rust_length']}, WDL={result['consensus']['wdl_length']}\n")
                f.write("\n")

                f.write("=== Stats Comparison ===\n")
                if 'error' in result['stats']:
                    f.write(f"Error: {result['stats']['error']}\n")
                else:
                    f.write("Comparing key statistics:\n")
                    f.write(f"Fields compared: {result['stats']['fields_compared']}\n")
                    f.write(
                        f"Fields within 1% tolerance or small absolute difference: {result['stats']['within_1_percent']}\n")
                    if result['stats']['within_1_percent_details']:
                        f.write("Fields within 1% tolerance or small absolute difference:\n")
                        for field, rust_val, wdl_val, diff in result['stats']['within_1_percent_details']:
                            if isinstance(diff, float):
                                f.write(f"  {field}: Rust={rust_val}, WDL={wdl_val}, Relative Diff={diff:.2f}%\n")
                            else:
                                f.write(f"  {field}: Rust={rust_val}, WDL={wdl_val}, Absolute Diff={diff}\n")
                    f.write(f"Fields within 5% tolerance (but not 1%): {result['stats']['within_5_percent']}\n")
                    if result['stats']['within_5_percent_details']:
                        f.write("Fields within 5% tolerance (but not 1%):\n")
                        for field, rust_val, wdl_val, diff in result['stats']['within_5_percent_details']:
                            f.write(f"  {field}: Rust={rust_val}, WDL={wdl_val}, Relative Diff={diff:.2f}%\n")
                    f.write(
                        f"Fields outside 5% tolerance or non-numeric mismatches: {result['stats']['outside_5_percent']}\n")
                    if result['stats']['outside_5_percent_details']:
                        f.write("Fields outside 5% tolerance or non-numeric mismatches:\n")
                        for field, rust_val, wdl_val, diff in result['stats']['outside_5_percent_details']:
                            if diff is None:
                                f.write(f"  {field}: Rust={rust_val}, WDL={wdl_val}\n")
                            else:
                                f.write(f"  {field}: Rust={rust_val}, WDL={wdl_val}, Relative Diff={diff:.2f}%\n")
                    f.write(f"Allele count mismatches: {result['stats']['allele_mismatches']}\n")
                    if result['stats']['allele_mismatch_details']:
                        f.write("Allele count mismatches:\n")
                        for key, rust_count, wdl_count, diff in result['stats']['allele_mismatch_details']:
                            f.write(f"  {key}: Rust={rust_count}, WDL={wdl_count}, Absolute Diff={diff}\n")
                    f.write("Note: 'coverage' array comparison skipped (requires specific handling)\n")
                    f.write("Note: WDL quantile keys were mapped from 'depth_q.XX' to 'depth_qXX' for comparison\n")
                    f.write(
                        "Note: Rust alignment metrics (mapped_reads, ercc_mapped_reads, ercc_mapped_paired) doubled for comparison due to paired-end concatenation\n")
                    f.write(
                        "Note: Paired-end metrics (mapped_paired, paired_inward, paired_outward, paired_other_orientation) excluded due to Rust pipeline concatenation\n")
                f.write("\n")

                f.write("=== Depth Plots Comparison ===\n")
                if result['depth_plots']['identical']:
                    f.write("Depth plots have identical file sizes\n")
                else:
                    f.write(
                        f"Depth plots differ in file size: Rust={result['depth_plots']['rust_size']}, WDL={result['depth_plots']['wdl_size']}\n")
                f.write("\n")

                f.write("=== QUAST Reports Comparison ===\n")
                f.write(f"Common keys: {result['quast_reports']['common_keys']}\n")
                f.write(f"Mismatches: {result['quast_reports']['mismatches']}\n")
                if result['quast_reports']['mismatch_details']:
                    f.write("Mismatched QUAST metrics:\n")
                    for key, rust_val, wdl_val in result['quast_reports']['mismatch_details']:
                        f.write(f"  {key}: Rust={rust_val}, WDL={wdl_val}\n")
                f.write("\n")

                f.write("=== Kraken Reports Comparison ===\n")
                f.write(f"Common keys: {result['kraken_reports']['common_keys']}\n")
                f.write(f"Mismatches: {result['kraken_reports']['mismatches']}\n")
                if result['kraken_reports']['mismatch_details']:
                    f.write("Mismatched Kraken metrics:\n")
                    for key, rust_val, wdl_val in result['kraken_reports']['mismatch_details']:
                        f.write(f"  {key}: Rust={rust_val}, WDL={wdl_val}\n")
                f.write("\n")

                f.write("=== ERCC Stats Comparison ===\n")
                f.write(f"Common keys: {result['ercc_stats']['common_keys']}\n")
                f.write(f"Mismatches: {result['ercc_stats']['mismatches']}\n")
                if result['ercc_stats']['mismatch_details']:
                    f.write("Mismatched ERCC metrics:\n")
                    for key, rust_val, wdl_val in result['ercc_stats']['mismatch_details']:
                        f.write(f"  {key}: Rust={rust_val}, WDL={wdl_val}\n")

            all_results.append(result)
        except Exception as e:
            print(f"Error processing {sample_name}: {str(e)}")

    with open("all_validation.tsv", 'w') as f:
        headers = [
            'Sample_name',
            'Variants_common_positions', 'Variants_rust_only_positions', 'Variants_wdl_only_positions',
            'Variants_mismatches',
            'Variants_mismatch_details', 'Variants_rust_only_details', 'Variants_wdl_only_details',
            'Consensus_rust_length', 'Consensus_wdl_length', 'Consensus_identical', 'Consensus_mismatch_position',
            'Consensus_mismatch_details',
            'Stats_fields_compared', 'Stats_within_1_percent', 'Stats_within_5_percent', 'Stats_outside_5_percent',
            'Stats_allele_mismatches',
            'Stats_within_1_percent_details', 'Stats_within_5_percent_details', 'Stats_outside_5_percent_details',
            'Stats_allele_mismatch_details',
            'Depth_plots_identical', 'Depth_plots_rust_size', 'Depth_plots_wdl_size',
            'QUAST_common_keys', 'QUAST_mismatches', 'QUAST_mismatch_details',
            'Kraken_common_keys', 'Kraken_mismatches', 'Kraken_mismatch_details',
            'ERCC_common_keys', 'ERCC_mismatches', 'ERCC_mismatch_details'
        ]
        f.write('\t'.join(headers) + '\n')

        for result in all_results:
            row = [result['sample_name']]
            if 'error' in result['variants']:
                row.extend(['error', 'error', 'error', 'error', 'error', 'error', 'error'])
            else:
                row.extend([
                    str(result['variants']['common_positions']),
                    str(result['variants']['rust_only_positions']),
                    str(result['variants']['wdl_only_positions']),
                    str(result['variants']['mismatches']),
                    ';'.join(
                        f"{key}:{','.join(f'REF={ref},ALT={','.join(alts)}' for ref, alts in rust_vars)}|{','.join(f'REF={ref},ALT={','.join(alts)}' for ref, alts in wdl_vars)}"
                        for key, rust_vars, wdl_vars in result['variants']['mismatch_details']),
                    ';'.join(
                        f"{key}:{','.join(f'REF={ref},ALT={','.join(alts)}' for ref, alts in vars)}" for key, vars in
                        result['variants']['rust_only_details']),
                    ';'.join(
                        f"{key}:{','.join(f'REF={ref},ALT={','.join(alts)}' for ref, alts in vars)}" for key, vars in
                        result['variants']['wdl_only_details'])
                ])
            if 'error' in result['consensus']:
                row.extend(['error', 'error', 'error', 'error', 'error'])
            else:
                row.extend([
                    str(result['consensus']['rust_length']),
                    str(result['consensus']['wdl_length']),
                    str(result['consensus']['identical']),
                    str(result['consensus']['mismatch_position'] or ''),
                    ';'.join(f"Pos={pos},Rust={rust},WDL={wdl}" for pos, rust, wdl in
                             result['consensus']['mismatch_details'])
                ])
            if 'error' in result['stats']:
                row.extend(['error', 'error', 'error', 'error', 'error', 'error', 'error', 'error', 'error'])
            else:
                row.extend([
                    str(result['stats']['fields_compared']),
                    str(result['stats']['within_1_percent']),
                    str(result['stats']['within_5_percent']),
                    str(result['stats']['outside_5_percent']),
                    str(result['stats']['allele_mismatches']),
                    ';'.join(
                        f"{field}:Rust={rust_val},WDL={wdl_val},Diff={diff if diff is None else (f'{diff:.2f}%' if isinstance(diff, float) else diff)}"
                        for field, rust_val, wdl_val, diff in result['stats']['within_1_percent_details']),
                    ';'.join(
                        f"{field}:Rust={rust_val},WDL={wdl_val},Diff={diff:.2f}%" for field, rust_val, wdl_val, diff in
                        result['stats']['within_5_percent_details']),
                    ';'.join(
                        f"{field}:Rust={rust_val},WDL={wdl_val},Diff={diff if diff is None else f'{diff:.2f}%'}" for
                        field, rust_val, wdl_val, diff in result['stats']['outside_5_percent_details']),
                    ';'.join(
                        f"{key}:Rust={rust_count},WDL={wdl_count},Diff={diff}" for key, rust_count, wdl_count, diff in
                        result['stats']['allele_mismatch_details'])
                ])
            row.extend([
                str(result['depth_plots']['identical']),
                str(result['depth_plots']['rust_size'] or ''),
                str(result['depth_plots']['wdl_size'] or '')
            ])
            row.extend([
                str(result['quast_reports']['common_keys']),
                str(result['quast_reports']['mismatches']),
                ';'.join(f"{key}:Rust={rust_val},WDL={wdl_val}" for key, rust_val, wdl_val in
                         result['quast_reports']['mismatch_details'])
            ])
            row.extend([
                str(result['kraken_reports']['common_keys']),
                str(result['kraken_reports']['mismatches']),
                ';'.join(f"{key}:Rust={rust_val},WDL={wdl_val}" for key, rust_val, wdl_val in
                         result['kraken_reports']['mismatch_details'])
            ])
            row.extend([
                str(result['ercc_stats']['common_keys']),
                str(result['ercc_stats']['mismatches']),
                ';'.join(f"{key}:Rust={rust_val},WDL={wdl_val}" for key, rust_val, wdl_val in
                         result['ercc_stats']['mismatch_details'])
            ])
            f.write('\t'.join(row) + '\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Validate multiple Rust pipeline results against WDL pipeline")
    parser.add_argument('input_file',
                        help='Whitespace-delimited file with columns <rust results dir> <wdl results dir>')

    args = parser.parse_args()
    validate_pipeline(args.input_file)