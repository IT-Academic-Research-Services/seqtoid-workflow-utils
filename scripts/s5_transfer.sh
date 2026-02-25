#!/usr/bin/env bash
# s5_transfer_low.sh - LOW concurrency test to avoid throttling

set -euo pipefail

INPUT_FILE="czid_s3_file_list.txt"
DEST_BUCKET="seqtoid-public-references"
DEST_PREFIX="phase1/"
NUMWORKERS=8             # LOW: 8 concurrent operations
CONCURRENCY=4            # LOW: 4 multipart parts per file
PART_SIZE=100            # Smaller parts to reduce per-request load
RETRY_COUNT=20           # More retries for flaky public access
LOG_FILE="czid_s5cmd_low_$(date +%Y%m%d_%H%M%S).log"
FAILED_FILE="czid_failed_low.txt"
CMD_FILE="czid_s5cmd_low_commands_$(date +%Y%m%d_%H%M%S).txt"

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

echo "=== CZ ID s5cmd LOW-CONCURRENCY TEST (throttling check) ==="
echo "Input file   : $INPUT_FILE ($(wc -l < "$INPUT_FILE") raw lines)"
echo "Dest         : s3://$DEST_BUCKET/$DEST_PREFIX"
echo "Workers      : $NUMWORKERS (very low)"
echo "Concurrency  : $CONCURRENCY (very low)"
echo "Part size    : $PART_SIZE MiB"
echo "Retries      : $RETRY_COUNT"
echo "Log          : $LOG_FILE"
echo "Failed list  : $FAILED_FILE"
echo "Dry run?     : $DRY_RUN"
echo "WARNING: This will be SLOW (~1-10 MB/s) but should avoid 403 throttling."
echo "================================================"

> "$CMD_FILE"

COUNT=0
while IFS= read -r line; do
    entry="${line#"${line%%[![:space:]]*}"}"   # ltrim
    entry="${entry%"${entry##*[![:space:]]}"}" # rtrim

    [[ -z "$entry" ]] && continue

    COUNT=$((COUNT + 1))
    entry="${entry%/}"

    src="s3://czid-public-references/$entry"
    dest="s3://$DEST_BUCKET/$DEST_PREFIX$entry"

    if [[ "$entry" == *.* || "$entry" == *.tar.gz || "$entry" == *.db || "$entry" == *.fa || "$entry" == *.fasta || "$entry" == *.json || "$entry" == *.txt || "$entry" == *.bed ]]; then
        src_pattern="$src"
        dest_path="$dest"
    else
        src_pattern="$src/*"
        dest_path="$dest/"
    fi

    echo "cp --no-clobber --concurrency $CONCURRENCY --part-size $PART_SIZE \"$src_pattern\" \"$dest_path\"" >> "$CMD_FILE"

    echo "Queued [$COUNT]: $entry → $src_pattern"
done < "$INPUT_FILE"

echo "Generated $COUNT commands → $CMD_FILE"

echo "First 5 commands in file:"
head -n 5 "$CMD_FILE" || echo "(empty?)"

if $DRY_RUN; then
    echo "Dry run complete."
    rm -f "$CMD_FILE"
    exit 0
fi

echo "Starting LOW-SPEED transfer (expect hours for 17 TB)..."
echo "Monitor: tail -f $LOG_FILE"
echo "Network: sudo apt install nload; nload eth0"

time s5cmd \
    --numworkers "$NUMWORKERS" \
    --retry-count "$RETRY_COUNT" \
    --log=debug \
    --stat \
    --no-sign-request \
    run "$CMD_FILE" 2>&1 | tee "$LOG_FILE"

echo "Transfer done. Exit code: $?"

grep -i -E '(ERROR|\+ERROR|failed|error:|NoSuchKey|AccessDenied|timeout|forbidden|retry|throttle)' "$LOG_FILE" > "$FAILED_FILE" || true

if [[ -s "$FAILED_FILE" ]]; then
    echo "Failures detected! See $FAILED_FILE"
    head -n 10 "$FAILED_FILE"
    echo "... (total: $(wc -l < "$FAILED_FILE"))"
else
    echo "No obvious failures in log — throttling likely avoided."
fi

echo "Clean up: rm $CMD_FILE"