#!/usr/bin/env bash
# aws_transfer_local.sh - Download anonymous, upload with creds via local temp

set -euo pipefail

INPUT_FILE="czid_s3_file_list.txt"
DEST_BUCKET="seqtoid-public-references"
DEST_PREFIX="phase1/"
PARALLEL_JOBS=8
TEMP_DIR="/tmp/s3_transfer"
LOG_FILE="czid_aws_local_transfer_$(date +%Y%m%d_%H%M%S).log"
FAILED_FILE="czid_aws_local_failed.txt"

DRY_RUN=false
while [[ $# -gt 0 ]]; do
    case $1 in
        --dry-run) DRY_RUN=true; shift ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done

DEST_PREFIX="${DEST_PREFIX%/}/"

if [[ ! -f "$INPUT_FILE" ]]; then
    echo "Error: Input file '$INPUT_FILE' not found"
    exit 1
fi

mkdir -p "$TEMP_DIR"

echo "=== CZ ID aws s3 cp via local temp (anonymous read, auth write) ==="
echo "Input file   : $INPUT_FILE ($(wc -l < "$INPUT_FILE") entries)"
echo "Dest         : s3://$DEST_BUCKET/$DEST_PREFIX"
echo "Parallel jobs: $PARALLEL_JOBS"
echo "Temp dir     : $TEMP_DIR"
echo "Log          : $LOG_FILE"
echo "Failed list  : $FAILED_FILE"
echo "Dry run?     : $DRY_RUN"
echo "WARNING: Will use local disk space for temp files."
echo "=================================================="

# Function to copy one entry
copy_entry() {
    entry="$1"
    entry="${entry%/}"

    src="s3://czid-public-references/$entry"
    dest="s3://$DEST_BUCKET/$DEST_PREFIX$entry"

    if [[ "$entry" == *.* || "$entry" == *.tar.gz || "$entry" == *.db || "$entry" == *.fa || "$entry" == *.fasta || "$entry" == *.json || "$entry" == *.txt || "$entry" == *.bed ]]; then
        src_path="$src"
        dest_path="$dest"
        is_dir=false
    else
        src_path="$src/"
        dest_path="$dest/"
        is_dir=true
    fi

    if $DRY_RUN; then
        echo "Would copy: $src_path to $dest_path (dir: $is_dir)"
        return 0
    fi

    local_temp="$TEMP_DIR/${entry##*/}"

    if $is_dir; then
        aws s3 cp --no-sign-request --recursive "$src_path" "$local_temp/" 2>&1
        if [[ $? -ne 0 ]]; then
            echo "Download failed: $entry" >> "$FAILED_FILE"
            return 1
        fi
        aws s3 cp --recursive "$local_temp/" "$dest_path" 2>&1
        if [[ $? -ne 0 ]]; then
            echo "Upload failed: $entry" >> "$FAILED_FILE"
            return 1
        fi
        rm -rf "$local_temp/"
    else
        aws s3 cp --no-sign-request "$src_path" "$local_temp" 2>&1
        if [[ $? -ne 0 ]]; then
            echo "Download failed: $entry" >> "$FAILED_FILE"
            return 1
        fi
        aws s3 cp "$local_temp" "$dest_path" 2>&1
        if [[ $? -ne 0 ]]; then
            echo "Upload failed: $entry" >> "$FAILED_FILE"
            return 1
        fi
        rm -f "$local_temp"
    fi
}

export -f copy_entry
export DRY_RUN TEMP_DIR DEST_BUCKET DEST_PREFIX FAILED_FILE

if $DRY_RUN; then
    cat "$INPUT_FILE" | parallel -j $PARALLEL_JOBS copy_entry {}
else
    cat "$INPUT_FILE" | parallel -j $PARALLEL_JOBS copy_entry {} 2>&1 | tee "$LOG_FILE"
fi

if [[ -s "$FAILED_FILE" ]]; then
    echo "Failures detected! See $FAILED_FILE"
    head -n 10 "$FAILED_FILE"
    echo "... (total: $(wc -l < "$FAILED_FILE"))"
else
    echo "No failures logged."
fi