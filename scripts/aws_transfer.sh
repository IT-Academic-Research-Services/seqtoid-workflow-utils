#!/usr/bin/env bash
# aws_transfer_local_v2.sh - Reliable cleanup + 403 logging

set -euo pipefail

INPUT_FILE="czid_s3_file_list.txt"
DEST_BUCKET="seqtoid-public-references"
DEST_PREFIX="phase1/"
PARALLEL_JOBS=4          # Lower to reduce disk pressure
LOG_FILE="czid_aws_local_v2_$(date +%Y%m%d_%H%M%S).log"
FAILED_FILE="czid_aws_local_failed_v2.txt"
TEMP_BASE="/tmp/s3_transfer_$$"  # Unique per run

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

mkdir -p "$TEMP_BASE"

echo "=== CZ ID aws s3 cp via local temp v2 (reliable cleanup) ==="
echo "Input file   : $INPUT_FILE ($(wc -l < "$INPUT_FILE") entries)"
echo "Dest         : s3://$DEST_BUCKET/$DEST_PREFIX"
echo "Parallel jobs: $PARALLEL_JOBS"
echo "Temp base    : $TEMP_BASE"
echo "Log          : $LOG_FILE"
echo "Failed list  : $FAILED_FILE"
echo "Dry run?     : $DRY_RUN"
echo "WARNING: Temp files cleaned per entry (even on failure)."
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
        echo "Would copy: $src_path → $dest_path (dir: $is_dir)"
        return 0
    fi

    # Per-entry temp dir (cleaned at end)
    local_temp=$(mktemp -d -p "$TEMP_BASE" "$entry.XXXXXXXX")
    trap 'rm -rf "$local_temp"' EXIT  # Cleanup on any exit

    # Check if dest already exists (skip download/upload)
    if aws s3 ls "$dest_path" >/dev/null 2>&1; then
        echo "Already exists: $dest_path" >> "$LOG_FILE"
        rm -rf "$local_temp"
        return 0
    fi

    if $is_dir; then
        aws s3 cp --no-sign-request --recursive "$src_path" "$local_temp/" 2>&1 || {
            echo "Download failed (403?): $entry" >> "$FAILED_FILE"
            rm -rf "$local_temp"
            return 1
        }
        aws s3 cp --recursive "$local_temp/" "$dest_path" 2>&1 || {
            echo "Upload failed: $entry" >> "$FAILED_FILE"
            rm -rf "$local_temp"
            return 1
        }
    else
        aws s3 cp --no-sign-request "$src_path" "$local_temp/file" 2>&1 || {
            if grep -q "AccessDenied\|403" <<< "$?"; then
                echo "403 on source: $entry" >> "$FAILED_FILE"
            else
                echo "Download failed: $entry" >> "$FAILED_FILE"
            fi
            rm -rf "$local_temp"
            return 1
        }
        aws s3 cp "$local_temp/file" "$dest_path" 2>&1 || {
            echo "Upload failed: $entry" >> "$FAILED_FILE"
            rm -rf "$local_temp"
            return 1
        }
    fi

    # Success — clean
    rm -rf "$local_temp"
    echo "Success: $entry" >> "$LOG_FILE"
}

export -f copy_entry
export DRY_RUN TEMP_BASE DEST_BUCKET DEST_PREFIX FAILED_FILE LOG_FILE

if $DRY_RUN; then
    cat "$INPUT_FILE" | parallel -j $PARALLEL_JOBS copy_entry {}
else
    cat "$INPUT_FILE" | parallel -j $PARALLEL_JOBS copy_entry {} 2>&1 | tee "$LOG_FILE"
fi

if [[ -s "$FAILED_FILE" ]]; then
    echo "Failures detected! See $FAILED_FILE"
    echo "403s for CZI escalation:"
    grep -i "403\|AccessDenied" "$FAILED_FILE" || echo "No 403s logged."
    echo "Total failed: $(wc -l < "$FAILED_FILE")"
else
    echo "No failures logged — all processed."
fi

# Final cleanup (in case of early exit)
rm -rf "$TEMP_BASE"
echo "Temp base cleaned: $TEMP_BASE"