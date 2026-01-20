import pandas as pd
import os
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

    # Placeholder for future steps
    # print("\n=== Step 3: Next Output Comparison (e.g., alignments) ===")
    # compare_alignments()

    # Add more steps as needed...

    print("\nComparison complete.")

if __name__ == '__main__':
    main()