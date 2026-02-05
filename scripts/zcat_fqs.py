import os
import subprocess
import sys
import argparse
from pathlib import Path




def main():
    parser = argparse.ArgumentParser(description="Compress .fq/.fastq files with pigz and delete originals")
    parser.add_argument("directory", nargs="?", default=".", help="Directory to scan (default: current)")
    parser.add_argument("--dry-run", action="store_true", help="Show what would happen, no changes")
    args = parser.parse_args()

    print("wat")
    root = Path(args.directory).resolve()
    if not root.is_dir():
        print(f"Error: {root} is not a directory")
        sys.exit(1)

    print(f"Scanning: {root}")
    if args.dry_run:
        print("DRY RUN — no files will be modified\n")

    extensions = ['.fq.gz', '.fastq.gz']

    file_pairs = []

    fqs = []
    for file in os.listdir(root):
        if file.endswith(extensions):
            fqs.append(file)
            print(f"using {file}")

    for file in fqs:
        file_elems = file.split(".")
        file_name = file_elems[0]
        name_cols = file_name.split("_")
        print(f"name cols {name_cols}")

        # this is dum and i dont care its staurday and im in a hirry
        pair_found = False
        for other_file in fqs:
            if pair_found == True:
                break
            other_file_elems = other_file.split(".")
            other_file_name = other_file_elems[0]
            other_name_cols = other_file_name.split("_")

            if name_cols[0] != other_name_cols[0]:
                if name_cols[1] == other_name_cols[1]:
                    if name_cols[3] == other_name_cols[3]:
                        if name_cols[2] == "500k" & other_name_cols[2] == "500k":
                            file_pairs.append((file, other_file))
                            pair_found = True
                            break
                        elif name_cols[2] == "900k" & other_name_cols[2] == "100k":
                            file_pairs.append((file, other_file))
                            pair_found = True
                            break

    for file_pair in file_pairs:
        print(f"{file_pair[0]} {file_pair[1]}")








