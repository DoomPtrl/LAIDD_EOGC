#!/usr/bin/env bash
set -euo pipefail

# Fix encoding
export LANG=C.UTF-8
export LC_ALL=C.UTF-8

###############################################################################
# CONFIGURATION
###############################################################################
RAW_ROOT="/data/Storage_DAS03/mentee/05_EOGC/proteomics/PDC000214/1/Raw Mass Spectra"
MZML_ROOT="/data/Storage_DAS03/home3/menteeF/yeongyu/Proteomics/01.mzml"
PARALLEL_JOBS=${PARALLEL_JOBS:-20}  # Default 20 for your 48-core system

###############################################################################
# CONVERT FUNCTION
###############################################################################
convert_one() {
  local RAW="$1"
  local rel=${RAW#"$RAW_ROOT"/}
  local sample_dir=$(dirname "$rel")
  local out_dir="$MZML_ROOT/$sample_dir/individual"
  local base=$(basename "$RAW" .raw)
  local out_file="$out_dir/${base}.mzML"
  
  # Skip if already exists
  if [[ -f "$out_file" ]]; then
    echo "[SKIP] $rel - already exists"
    return 0
  fi
  
  echo "[START] $rel"
  mkdir -p "$out_dir"
  
  # Convert
  ThermoRawFileParser \
    -i="$RAW" \
    -o="$out_dir" \
    -f=1
  
  if [[ -f "$out_file" ]]; then
    echo "[DONE] $rel"
  else
    echo "[FAIL] $rel" >> "$MZML_ROOT/failed.txt"
  fi
}

export -f convert_one
export RAW_ROOT MZML_ROOT

###############################################################################
# RUN
###############################################################################
echo "Converting $(find "$RAW_ROOT" -name '*.raw' | wc -l) RAW files with $PARALLEL_JOBS jobs"

if command -v parallel &>/dev/null; then
  find "$RAW_ROOT" -name '*.raw' -print0 | \
  parallel -0 -j "$PARALLEL_JOBS" convert_one
else
  find "$RAW_ROOT" -name '*.raw' | while read -r RAW; do
    convert_one "$RAW"
  done
fi

echo "Complete. Check $MZML_ROOT/failed.txt for any failures."