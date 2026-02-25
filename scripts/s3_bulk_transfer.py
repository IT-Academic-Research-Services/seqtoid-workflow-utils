#!/usr/bin/env python3
"""
s3_bulk_copy_boto3.py
Copy objects from czid-public-references → seqtoid-public-references/phase1/
Preserves full key structure.
Server-side copy → fast & cheap.
Parallelized with ThreadPoolExecutor.
Run with: python3 s3_bulk_copy_boto3.py
"""

import boto3
from concurrent.futures import ThreadPoolExecutor, as_completed
from botocore.exceptions import ClientError
import sys
from datetime import datetime

# ────────────────────────────────────────────────
# CONFIG
# ────────────────────────────────────────────────
SOURCE_BUCKET = "czid-public-references"
DEST_BUCKET   = "seqtoid-public-references"
DEST_PREFIX   = "phase1/"   # will be prepended to every source key

# List of prefixes/keys to copy (your original list, one per line)
# If you prefer, read from file instead of hardcoding
PREFIXES_TO_COPY = [
    "host_filter",
    "taxonomy",
    "ncbi-indexes-prod/2022-06-02/index-generation-2/accession2taxid.marisa",
    # ... add ALL 78 entries here (copy-paste from your s5out.txt "Queued [n]: ..." lines)
    # For brevity I'm showing first few — replace with your full list!
    "test/viral-alignment-indexes/viral_nt_loc.marisa",
    # etc.
]

# Or better: read from your original input file
INPUT_FILE = "czid_s3_file_list.txt"  # same as your s5cmd script

MAX_WORKERS = 100  # tune: 50–300 depending on network/CPU; 100 is safe start
# ────────────────────────────────────────────────

def load_prefixes():
    if INPUT_FILE and INPUT_FILE.strip():
        try:
            with open(INPUT_FILE, "r") as f:
                lines = [line.strip() for line in f if line.strip()]
            print(f"Loaded {len(lines)} prefixes from {INPUT_FILE}")
            return lines
        except FileNotFoundError:
            print(f"File {INPUT_FILE} not found → falling back to hardcoded list")
    return PREFIXES_TO_COPY

def copy_object(s3_client, src_key, dest_key):
    try:
        # Server-side copy
        s3_client.copy_object(
            Bucket=DEST_BUCKET,
            CopySource={'Bucket': SOURCE_BUCKET, 'Key': src_key},
            Key=dest_key,
            # Optional: add MetadataDirective='REPLACE' if you want to change metadata
            # ACL='private',  # default is private anyway
        )
        print(f"OK  {src_key} → {dest_key}")
        return True, None
    except ClientError as e:
        err = e.response['Error']
        print(f"ERR {src_key} → {err['Code']} {err['Message']}")
        return False, src_key

def main():
    s3 = boto3.client('s3')  # assumes IAM role / credentials are set

    prefixes = load_prefixes()
    if not prefixes:
        print("No prefixes to copy. Exiting.")
        sys.exit(1)

    print(f"\nStarting bulk server-side copy")
    print(f"  Source bucket: {SOURCE_BUCKET}")
    print(f"  Dest   bucket: {DEST_BUCKET}/{DEST_PREFIX}")
    print(f"  Items: {len(prefixes)}")
    print(f"  Parallel workers: {MAX_WORKERS}")
    print("-" * 60)

    failed = []
    total = 0
    copied = 0
    start = datetime.now()

    # If prefix ends with / → it's a "directory" → copy recursively
    # Otherwise → single object
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = []

        for prefix in prefixes:
            prefix = prefix.rstrip('/')  # clean

            if prefix.endswith('/'):  # treat as dir prefix
                # List objects under this prefix
                paginator = s3.get_paginator('list_objects_v2')
                for page in paginator.paginate(Bucket=SOURCE_BUCKET, Prefix=prefix):
                    if 'Contents' not in page:
                        continue
                    for obj in page['Contents']:
                        src_key = obj['Key']
                        if src_key.endswith('/'):  # skip "folder" markers
                            continue
                        dest_key = DEST_PREFIX + src_key
                        futures.append(executor.submit(copy_object, s3, src_key, dest_key))
                        total += 1
            else:
                # single file
                dest_key = DEST_PREFIX + prefix
                futures.append(executor.submit(copy_object, s3, prefix, dest_key))
                total += 1

        # Wait & collect results
        for future in as_completed(futures):
            success, failed_key = future.result()
            if success:
                copied += 1
            else:
                if failed_key:
                    failed.append(failed_key)

    duration = datetime.now() - start
    print("\n" + "="*60)
    print(f"Finished in {duration}")
    print(f"Copied:  {copied}/{total}")
    print(f"Failed:  {len(failed)}")

    if failed:
        print("\nFailed keys:")
        for k in failed[:20]:
            print(k)
        if len(failed) > 20:
            print(f"... +{len(failed)-20} more")
        with open("failed_copies.txt", "w") as f:
            f.write("\n".join(failed))
        print("→ Saved to failed_copies.txt")

if __name__ == "__main__":
    main()