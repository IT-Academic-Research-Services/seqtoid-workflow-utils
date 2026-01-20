import pandas as pd
import os
import glob
from typing import List, Dict

# Define paths to the directories
CZID_DIR = 'czid'
SEQTOID_DIR = 'seqtoid'

# List of expected sample IDs (from the provided NCBI IDs)
EXPECTED_SAMPLES = [
    'ERR11417004',
    'SRR23038836',
    'SRR11454628',
    'SRR12876565',
    'SRR14579537',
    'SRR10903401',
    'SRR1304850',
    'SRR13227005',
    'SRR13227004',
    'SRR13227003',
    'SRR11278904'
]

def compare_metadata():
    """
    Step 1: Compare sample_metadata.csv files from czid and seqtoid directories.
    - Loads the CSVs using pandas.
    - Sorts by 'sample_name' to handle potential order differences.
    - Checks for exact equality.
    - If not equal, computes and prints differences.
    - Also checks for missing/extra samples relative to the expected list.
    """
    metadata_file = 'sample_metadata.csv'

    czid_path = os.path.join(CZID_DIR, metadata_file)
    seqtoid_path = os.path.join(SEQTOID_DIR, metadata_file)

    if not os.path.exists(czid_path):
        print(f"Error: {czid_path} not found.")
        return
    if not os.path.exists(seqtoid_path):
        print(f"Error: {seqtoid_path} not found.")
        return

    # Read CSVs
    czid_df = pd.read_csv(czid_path)
    seqtoid_df = pd.read_csv(seqtoid_path)

    # Sort by sample_name for comparison
    czid_df_sorted = czid_df.sort_values('sample_name').reset_index(drop=True)
    seqtoid_df_sorted = seqtoid_df.sort_values('sample_name').reset_index(drop=True)

    print("Comparing sample_metadata.csv...")

    if czid_df_sorted.equals(seqtoid_df_sorted):
        print("  - Metadata files are identical.")
    else:
        print("  - Metadata files differ.")
        # Compute differences
        diff = czid_df_sorted.compare(seqtoid_df_sorted)
        if not diff.empty:
            print("    Differences (czid vs seqtoid):")
            print(diff)
        else:
            print("    Differences in structure (e.g., row count or columns).")
            if list(czid_df_sorted.columns) != list(seqtoid_df_sorted.columns):
                print("      Column mismatch:")
                print(f"        czid: {czid_df_sorted.columns}")
                print(f"        seqtoid: {seqtoid_df_sorted.columns}")
            if len(czid_df_sorted) != len(seqtoid_df_sorted):
                print(f"      Row count mismatch: czid={len(czid_df_sorted)}, seqtoid={len(seqtoid_df_sorted)}")

    # Check for missing/extra samples in each dataframe relative to expected list
    # Normalize expected samples (in case of suffixes, but compare base names if needed; for now, exact match)
    for df_name, df in [('czid', czid_df), ('seqtoid', seqtoid_df)]:
        actual_samples = set(df['sample_name'].astype(str))
        expected_set = set(EXPECTED_SAMPLES)

        missing = expected_set - actual_samples
        extra = actual_samples - expected_set

        if missing:
            print(f"  - Missing samples in {df_name}: {', '.join(sorted(missing))}")
        if extra:
            print(f"  - Extra samples in {df_name}: {', '.join(sorted(extra))}")

def compare_overviews():
    """
    Step 2: Compare sample_overviews.csv files from czid and seqtoid directories.
    - Loads the CSVs using pandas.
    - Sorts by 'sample_name' to handle potential order differences.
    - Checks for exact equality.
    - If not equal, computes and prints differences.
    - Also checks for missing/extra samples relative to the expected list.
    """
    overviews_file = 'sample_overviews.csv'

    czid_path = os.path.join(CZID_DIR, overviews_file)
    seqtoid_path = os.path.join(SEQTOID_DIR, overviews_file)

    if not os.path.exists(czid_path):
        print(f"Error: {czid_path} not found.")
        return
    if not os.path.exists(seqtoid_path):
        print(f"Error: {seqtoid_path} not found.")
        return

    # Read CSVs, allowing for mixed types (e.g., empty strings in numeric columns)
    czid_df = pd.read_csv(czid_path, dtype=str)  # Read as string to handle mixed types consistently
    seqtoid_df = pd.read_csv(seqtoid_path, dtype=str)

    # Sort by sample_name for comparison
    czid_df_sorted = czid_df.sort_values('sample_name').reset_index(drop=True)
    seqtoid_df_sorted = seqtoid_df.sort_values('sample_name').reset_index(drop=True)

    print("Comparing sample_overviews.csv...")

    if czid_df_sorted.equals(seqtoid_df_sorted):
        print("  - Overviews files are identical.")
    else:
        print("  - Overviews files differ.")
        # Compute differences
        diff = czid_df_sorted.compare(seqtoid_df_sorted)
        if not diff.empty:
            print("    Differences (czid vs seqtoid):")
            print(diff)
        else:
            print("    Differences in structure (e.g., row count or columns).")
            if list(czid_df_sorted.columns) != list(seqtoid_df_sorted.columns):
                print("      Column mismatch:")
                print(f"        czid: {czid_df_sorted.columns}")
                print(f"        seqtoid: {seqtoid_df_sorted.columns}")
            if len(czid_df_sorted) != len(seqtoid_df_sorted):
                print(f"      Row count mismatch: czid={len(czid_df_sorted)}, seqtoid={len(seqtoid_df_sorted)}")

    # Check for missing/extra samples in each dataframe relative to expected list
    # Normalize expected samples (in case of suffixes, but compare base names if needed; for now, exact match)
    for df_name, df in [('czid', czid_df), ('seqtoid', seqtoid_df)]:
        actual_samples = set(df['sample_name'].astype(str))
        expected_set = set(EXPECTED_SAMPLES)

        missing = expected_set - actual_samples
        extra = actual_samples - expected_set

        if missing:
            print(f"  - Missing samples in {df_name}: {', '.join(sorted(missing))}")
        if extra:
            print(f"  - Extra samples in {df_name}: {', '.join(sorted(extra))}")

def compare_taxon_reports():
    """
    Step 3: Compare per-sample taxon_report.csv files from czid and seqtoid directories.
    - For each expected sample, finds files matching {sample}_*_taxon_report.csv using glob.
    - Assumes exactly one file per sample per directory.
    - Loads the CSVs using pandas.
    - Sorts by 'tax_id' (which can be negative) to handle potential order differences.
    - Checks for exact equality per sample.
    - If not equal, computes and prints differences.
    - Reports missing files for samples.
    """
    print("Comparing per-sample taxon_report.csv files...")

    missing_in_czid = []
    missing_in_seqtoid = []

    for sample in EXPECTED_SAMPLES:
        # Find files in czid
        czid_pattern = os.path.join(CZID_DIR, f"{sample}_*_taxon_report.csv")
        czid_files = glob.glob(czid_pattern)
        if len(czid_files) != 1:
            print(f"  - For sample {sample} in czid: {'No file found' if len(czid_files) == 0 else 'Multiple files found'}")
            missing_in_czid.append(sample)
            continue
        czid_path = czid_files[0]

        # Find files in seqtoid
        seqtoid_pattern = os.path.join(SEQTOID_DIR, f"{sample}_*_taxon_report.csv")
        seqtoid_files = glob.glob(seqtoid_pattern)
        if len(seqtoid_files) != 1:
            print(f"  - For sample {sample} in seqtoid: {'No file found' if len(seqtoid_files) == 0 else 'Multiple files found'}")
            missing_in_seqtoid.append(sample)
            continue
        seqtoid_path = seqtoid_files[0]

        # Read CSVs, allowing for mixed types
        czid_df = pd.read_csv(czid_path, dtype=str)
        seqtoid_df = pd.read_csv(seqtoid_path, dtype=str)

        # Sort by tax_id for comparison (convert to int for proper sorting, including negatives)
        czid_df['tax_id'] = pd.to_numeric(czid_df['tax_id'], errors='coerce')
        seqtoid_df['tax_id'] = pd.to_numeric(seqtoid_df['tax_id'], errors='coerce')
        czid_df_sorted = czid_df.sort_values('tax_id').reset_index(drop=True)
        seqtoid_df_sorted = seqtoid_df.sort_values('tax_id').reset_index(drop=True)

        # Convert back to str for comparison if needed, but since dtype=str initially, and equals handles it
        print(f"  - Comparing for sample {sample}...")

        if czid_df_sorted.equals(seqtoid_df_sorted):
            print("    - Taxon report files are identical.")
        else:
            print("    - Taxon report files differ.")
            # Compute differences
            diff = czid_df_sorted.compare(seqtoid_df_sorted)
            if not diff.empty:
                print("      Differences (czid vs seqtoid):")
                print(diff)
            else:
                print("      Differences in structure (e.g., row count or columns).")
                if list(czid_df_sorted.columns) != list(seqtoid_df_sorted.columns):
                    print("        Column mismatch:")
                    print(f"          czid: {czid_df_sorted.columns}")
                    print(f"          seqtoid: {seqtoid_df_sorted.columns}")
                if len(czid_df_sorted) != len(seqtoid_df_sorted):
                    print(f"        Row count mismatch: czid={len(czid_df_sorted)}, seqtoid={len(seqtoid_df_sorted)}")

    if missing_in_czid:
        print(f"  - Missing taxon reports in czid for samples: {', '.join(missing_in_czid)}")
    if missing_in_seqtoid:
        print(f"  - Missing taxon reports in seqtoid for samples: {', '.join(missing_in_seqtoid)}")

def main():
    """
    Main entry point for the comparison script.
    This script is structured to allow adding more comparison steps in sequence.
    Each step compares a specific output artifact from the CZID pipeline.
    """
    print("Starting CZID pipeline output comparison...")

    print("\n=== Step 1: Sample Metadata Comparison ===")
    compare_metadata()

    print("\n=== Step 2: Sample Overviews Comparison ===")
    compare_overviews()

    print("\n=== Step 3: Sample Taxon Reports Comparison ===")
    compare_taxon_reports()

    # Placeholder for future steps
    # print("\n=== Step 4: Next Output Comparison ===")
    # compare_next()

    # Add more steps as needed...

    print("\nComparison complete.")

if __name__ == '__main__':
    main()