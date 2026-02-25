import os
import re
import logging
from typing import Set

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def extract_s3_keys(directory: str, bucket_name: str = 'czid-public-references') -> Set[str]:
    """
    Extract only concrete S3 object keys (must contain '/') from the repo.
    Skips:
    - .md files entirely
    - Full-line comments
    - Inline comments
    - Top-level prefixes (no '/' after bucket → likely just base dirs)
    """
    # Pattern requires at least one / after the bucket (real file/subpath)
    s3_pattern = re.compile(rf's3://{bucket_name}/([^ \t\n\r\f\v"\'\`]+/[^ \t\n\r\f\v"\'\`]+)')
    # ↑ this ensures we match something like bucket/dir/... or bucket/file (but not bucket alone)

    unique_keys: Set[str] = set()

    for root, _, files in os.walk(directory):
        for file_name in files:
            if file_name.lower().endswith('.md'):
                logger.debug(f"Skipping Markdown file: {os.path.join(root, file_name)}")
                continue

            file_path = os.path.join(root, file_name)
            try:
                with open(file_path, 'r', encoding='utf-8') as f:
                    for line_num, line in enumerate(f, 1):
                        stripped = line.strip()

                        # Skip full-line comments
                        if stripped.startswith('#'):
                            continue

                        # Strip inline comments
                        if '#' in line:
                            line = line.split('#', 1)[0]

                        matches = s3_pattern.findall(line)
                        for key in matches:
                            # Clean trailing junk (punctuation, quotes, extra slash)
                            key = key.rstrip('.,;:"\'/')
                            key = key.strip()
                            if key and '/' in key:  # double-check (redundant but safe)
                                unique_keys.add(key)
                                logger.info(f"Found real file key '{key}' in {file_path}:{line_num}")
            except UnicodeDecodeError:
                logger.debug(f"Skipping binary/non-UTF8 file: {file_path}")
            except Exception as e:
                logger.warning(f"Error reading {file_path}: {str(e)}")

    return unique_keys

def generate_file_list(output_file: str, keys: Set[str]):
    with open(output_file, 'w') as f:
        for key in sorted(keys):
            f.write(f"{key}\n")
    logger.info(f"Generated clean file list at {output_file} with {len(keys)} entries")

def main():
    repo_directory = '.'  # ← change if repo is not current dir
    output_file    = 's3_file_list_real_files.txt'

    keys = extract_s3_keys(repo_directory)
    generate_file_list(output_file, keys)

if __name__ == "__main__":
    main()