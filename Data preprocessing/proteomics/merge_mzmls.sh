#!/bin/bash
set -euo pipefail

MZML_ROOT="/data/Storage_DAS03/home3/menteeF/yeongyu/Proteomics/01.mzml"

echo "Starting mzML merging..."

# Exclude any path that contains a single quote
find "$MZML_ROOT" -type d -name individual ! -path "*'*" -print0 |
while IFS= read -r -d '' INDIV_DIR; do
  echo "Ì†ΩÌ¥ç Checking directory: $INDIV_DIR"

  # Skip empty dirs
  shopt -s nullglob
  files=( "$INDIV_DIR"/*.mzML )
  shopt -u nullglob
  if (( ${#files[@]} == 0 )); then
    echo "‚ö†Ô∏è  No mzML files in $INDIV_DIR, skipping."
    continue
  fi

  PROPRIETARY_DIR=$(dirname "$INDIV_DIR")
  SAMPLE_DIR=$(dirname "$PROPRIETARY_DIR")
  SAMPLE_NAME=$(basename "$SAMPLE_DIR")
  CORE_NAME=$(echo "$SAMPLE_NAME" | cut -d'_' -f2,3)
  OUT_FILE="$PROPRIETARY_DIR/${CORE_NAME}_merged.mzML"

  echo "Ì†ΩÌ≥¶ Merging: ${#files[@]} files ‚Üí $OUT_FILE"
  FileMerger -in "${files[@]}" -out "$OUT_FILE"
  echo "‚úÖ [$SAMPLE_NAME] merged to $OUT_FILE"
done

echo "‚úÖ All directories processed."
