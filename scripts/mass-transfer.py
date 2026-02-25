import boto3
import os
import logging
import concurrent.futures
from botocore.exceptions import ClientError
import tempfile
import csv
import argparse

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# S3 client
s3 = boto3.client('s3')

# Source and destination details
source_bucket = 'czid-public-references'
dest_bucket = 'seqtoid-public-references'
dest_prefix = 'phase1/'

def file_exists_in_s3(bucket, key):
    try:
        s3.head_object(Bucket=bucket, Key=key)
        return True
    except ClientError as e:
        if e.response['Error']['Code'] == '404':
            return False
        raise

def collect_concrete_keys(bucket, entry):
    keys = []
    # Check if exact key is a file
    try:
        s3.head_object(Bucket=bucket, Key=entry)
        keys.append(entry)
    except ClientError as e:
        if e.response['Error']['Code'] == '404':
            # Assume directory: list under prefix (add / if needed)
            prefix = entry if entry.endswith('/') else entry + '/'
            paginator = s3.get_paginator('list_objects_v2')
            for page in paginator.paginate(Bucket=bucket, Prefix=prefix):
                for obj in page.get('Contents', []):
                    if not obj['Key'].endswith('/'):  # Skip any placeholder dirs
                        keys.append(obj['Key'])
        else:
            raise
    return keys

def process_file(key, results):
    source_path = f"s3://{source_bucket}/{key}"
    dest_key = os.path.join(dest_prefix, key).replace('\\', '/')  # Normalize
    dest_path = f"s3://{dest_bucket}/{dest_key}"
    logger.info(f"Processing file: {key} (dest: {dest_key})")

    if file_exists_in_s3(dest_bucket, dest_key):
        logger.info(f"File already exists in destination: {dest_key}")
        results.append((source_path, dest_path, ""))
        return True

    local_path = None
    try:
        # Create a temporary file
        with tempfile.NamedTemporaryFile(delete=False) as temp_file:
            local_path = temp_file.name

        # Download from source
        s3.download_file(source_bucket, key, local_path)
        logger.info(f"Downloaded {key} to local: {local_path}")

        # Upload to destination
        s3.upload_file(local_path, dest_bucket, dest_key)
        logger.info(f"Uploaded {key} to {dest_key}")

        # Delete local file
        os.remove(local_path)
        logger.info(f"Deleted local file: {local_path}")

        logger.info(f"Successfully processed {key}")
        results.append((source_path, dest_path, ""))
        return True
    except Exception as e:
        error_msg = str(e)
        logger.error(f"Failed to process {key}: {error_msg}")
        if local_path and os.path.exists(local_path):
            os.remove(local_path)
        results.append((source_path, "FAILED", error_msg))
        return False

def main(input_file, output_file):
    # Read entries
    with open(input_file, 'r') as f:
        entries = [line.strip() for line in f if line.strip()]

    # Collect all concrete keys (resolving dirs)
    all_keys = set()
    for entry in entries:
        keys = collect_concrete_keys(source_bucket, entry)
        logger.info(f"Resolved {entry} to {len(keys)} concrete files")
        all_keys.update(keys)

    results = []
    successful = 0
    failed = 0

    with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:  # Adjust as needed
        futures = [executor.submit(process_file, key, results) for key in sorted(all_keys)]
        for future in concurrent.futures.as_completed(futures):
            if future.result():
                successful += 1
            else:
                failed += 1

    # Write results to CSV (sorted by source)
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['source s3 file', 'destination s3 file (or FAILED)', 'reason for failure'])
        writer.writerows(sorted(results))

    logger.info(f"Processing complete. Successful: {successful}, Failed: {failed}. Output: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Copy S3 files/dirs from list")
    parser.add_argument('--input_file', default='s3_file_list_final.txt', help='Path to input list')
    parser.add_argument('--output_file', default='transfer_results.csv', help='Path to output CSV')
    args = parser.parse_args()
    main(args.input_file, args.output_file)