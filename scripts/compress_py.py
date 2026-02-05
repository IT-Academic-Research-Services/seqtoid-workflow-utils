#!/usr/bin/env python3
"""
Compress all .fq / .fastq files in a directory using pigz and delete originals.

Usage:
    python compress_fq.py                  # current directory
    python compress_fq.py /path/to/dir     # specific directory
    python compress_fq.py --dry-run        # preview only, no changes
"""

import os
import subprocess
import sys
import argparse
from pathlib import Path


def compress_and_delete(fq_path: Path, dry_run: bool = False):
    gz_path = fq_path.with_suffix(fq_path.suffix + '.gz')

    if gz_path.exists():
        print(f"SKIP: {gz_path} already exists → skipping {fq_path}")
        return

    print(f"Compressing: {fq_path} → {gz_path}")

    if dry_run:
        return

    try:
        # Run pigz (multi-threaded gzip) and capture output
        result = subprocess.run(
            ['pigz', '-p', '8', '-c', str(fq_path)],   # -p 8 = 8 threads; adjust as needed
            check=True,
            stdout=open(gz_path, 'wb'),
            stderr=subprocess.PIPE,
            text=True
        )

        if result.returncode == 0:
            print(f"Success: {gz_path} created")
            # Delete original only after successful compression
            fq_path.unlink()
            print(f"Deleted original: {fq_path}")
        else:
            print(f"pigz failed for {fq_path}: {result.stderr}")
            # Optional: remove partial .gz if failed
            if gz_path.exists():
                gz_path.unlink()
                print(f"Removed partial: {gz_path}")

    except subprocess.CalledProcessError as e:
        print(f"Error compressing {fq_path}: {e}")
        if gz_path.exists():
            gz_path.unlink()
    except FileNotFoundError:
        print("pigz not found — install with: sudo apt install pigz")
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description="Compress .fq/.fastq files with pigz and delete originals")
    parser.add_argument("directory", nargs="?", default=".", help="Directory to scan (default: current)")
    parser.add_argument("--dry-run", action="store_true", help="Show what would happen, no changes")
    parser.add_argument("--include-fastq", action="store_true", help="Also process .fastq files")
    args = parser.parse_args()

    root = Path(args.directory).resolve()
    if not root.is_dir():
        print(f"Error: {root} is not a directory")
        sys.exit(1)

    print(f"Scanning: {root}")
    if args.dry_run:
        print("DRY RUN — no files will be modified\n")

    extensions = ['.fq']
    if args.include_fastq:
        extensions.append('.fastq')

    found = 0
    for ext in extensions:
        for fq in root.glob(f"*{ext}"):
            found += 1
            compress_and_delete(fq, dry_run=args.dry_run)

    if found == 0:
        print(f"No .{','.join(extensions)} files found in {root}")
    else:
        print(f"\nProcessed {found} files")


if __name__ == "__main__":
    main()