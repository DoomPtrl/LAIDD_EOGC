#!/bin/bash
#SBATCH -J minsea_raw2mzml
#SBATCH --nodelist=n2
#SBATCH --partition=node2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --output=slurm-%j.out

RAW_ROOT="./raw"
MZML_ROOT="./mzml"
THREADS=30

echo " Raw 2 mzML (Start)"

#i 모든 .raw 파일 찾기
find "$RAW_ROOT" -type f -name "*.raw" > raw_list.txt

# 병렬 변환
cat raw_list.txt | parallel -j $THREADS '
  RAW_FILE="{}"
  REL_PATH="${RAW_FILE#'$RAW_ROOT'/}"               # raw_root 기준 상대 경로
  RAW_NAME=$(basename "$RAW_FILE")
  TARGET_DIR="'$MZML_ROOT'/$(dirname "$REL_PATH")/individual"

  mkdir -p "$TARGET_DIR"
  echo "   🔸 Converting: $RAW_NAME → $TARGET_DIR"
  ThermoRawFileParser -i "$RAW_FILE" -o "$TARGET_DIR" -f 1
'


echo " Raw 2 mzML (Done)"