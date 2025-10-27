#!/usr/bin/env bash
# Split into ONE file per job and submit; resume-safe.
set -eu; (set -o 2>/dev/null | grep -q pipefail) && set -o pipefail

ORIG_JSON="JSON/global_250715.fixed.json"      # base config with full mzml_paths
PATH_TXT="JSON/mzmlPath.txt"                   # one CSV line of "path1","path2",...
BATCH_SIZE=1                                   # one file per job
OUTPUT_DIR="JSON/parts_single"                 # per-file JSONs
OUTPUT_SAGE_DIR="02.sage/Proteome_single"      # outputs per file
SAGE_BIN="/data/Storage_DAS03/home3/menteeF/.conda/envs/mzml_env/bin/sage"
PARTITION="node2"
CPUS=8                                         # threads per job (lower RAM); bump later if stable

mkdir -p "$OUTPUT_DIR" "$OUTPUT_SAGE_DIR"

# Make mzmlPath.txt if missing (CSV from ORIG_JSON)
if [[ ! -s "$PATH_TXT" ]]; then
  jq -r '.mzml_paths | @csv' "$ORIG_JSON" > "$PATH_TXT"
fi

# Read CSV -> array, strip quotes
line="$(tr -d '\r' < "$PATH_TXT" | head -n1)"
IFS=',' read -r -a mzml_paths <<< "$line"
for i in "${!mzml_paths[@]}"; do mzml_paths[$i]="${mzml_paths[$i]//\"/}"; done
TOTAL=${#mzml_paths[@]}
echo "[i] Found $TOTAL mzMLs"

count=1
for ((i=0; i<TOTAL; i+=BATCH_SIZE)); do
  part_paths=( "${mzml_paths[@]:i:BATCH_SIZE}" )
  jq_array=$(printf '%s\n' "${part_paths[@]}" | jq -R . | jq -s .)

  idx=$(printf "%03d" "$count")
  out_json="${OUTPUT_DIR}/single_${idx}.json"

  # Write per-file JSON + lighten memory: max_peaks=120, max_variable_mods=1
  jq --argjson paths "$jq_array" \
     '.mzml_paths = $paths | .max_peaks=120 | .database.max_variable_mods=1' \
     "$ORIG_JSON" > "$out_json"
  echo "[i] Wrote $out_json"

  # Resume-safe: skip if already finished
  if [[ -s "${OUTPUT_SAGE_DIR}/res_${idx}/results.sage.tsv" ]]; then
    echo "[i] Skipping ${idx} (already done)"
    ((count++)); continue
  fi

  sbatch <<EOT
#!/bin/bash
#SBATCH -J sage_single_${idx}
#SBATCH -p ${PARTITION}
#SBATCH -n 1
#SBATCH --cpus-per-task=${CPUS}
#SBATCH -t 1-12:00:00
#SBATCH --mem=0
#SBATCH -o .single_${idx}_%j.out
#SBATCH -e .single_${idx}_%j.err

set -euo pipefail
ulimit -s unlimited
export RUST_MIN_STACK=\$((32*1024*1024))
export RAYON_NUM_THREADS=${CPUS}
export OMP_NUM_THREADS=${CPUS}
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export MALLOC_ARENA_MAX=1

OUT="${OUTPUT_SAGE_DIR}/res_${idx}"
mkdir -p "\$OUT"
"${SAGE_BIN}" "${out_json}" -o "\$OUT"
EOT

  ((count++))
done

echo "[i] Submitted \$((count-1)) single-file jobs (CPUS=${CPUS})"
