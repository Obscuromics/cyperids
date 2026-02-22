#!/bin/bash
# Script to build round-1-4 cumulative libraries for RepeatMasker

set -euo pipefail

# Parent directory containing RM_* directories
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <parent_directory>"
    exit 1
fi

PARENT_DIR="$1"
OUT_DIR="$PARENT_DIR"

# Loop through all RM_* directories
for RM_DIR in "$PARENT_DIR"/RM_*; do
    [[ -d "$RM_DIR" ]] || continue
    echo "Processing $RM_DIR"

    LOG_FILE="$RM_DIR/rmod.log"
    if [[ ! -f "$LOG_FILE" ]]; then
        echo "No rmod.log found in $RM_DIR, skipping..."
        continue
    fi

    # Extract species name from rmod.log (after ID_)
    SPECIES=$(sed -n 's/.*ID_\([^[:space:]]*\).*/\1/p' "$LOG_FILE" | head -n1)

    if [[ -z "$SPECIES" ]]; then
        echo "Could not extract species name from $LOG_FILE, skipping..."
        continue
    fi

    echo "  Species identified as: $SPECIES"

    # Temp file to store concatenation
    TEMP_CONCAT=$(mktemp)

    # Loop over rounds 1–4
    for ROUND in round-{1..4}; do
        CONSENSI_FILE="$RM_DIR/$ROUND/consensi.fa"
        if [[ -s "$CONSENSI_FILE" ]]; then
            cat "$CONSENSI_FILE" >> "$TEMP_CONCAT"
            printf '\n' >> "$TEMP_CONCAT"
        else
            echo "Warning: $CONSENSI_FILE not found or empty, skipping"
        fi
    done

    # Output file
    OUT_FILE="$OUT_DIR/round-1-2-3-4-consensi_${SPECIES}.fa"

    # Remove duplicates with CD-HIT (95% identity)
    cd-hit-est -i "$TEMP_CONCAT" -o "$OUT_FILE" -c 0.95 -n 10

    echo "Saved deduplicated cumulative library to $OUT_FILE"

    rm "$TEMP_CONCAT"
done
