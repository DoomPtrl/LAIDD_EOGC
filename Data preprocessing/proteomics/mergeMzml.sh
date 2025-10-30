#!/bin/bash
#SBATCH -J minsea_raw2mzml
#SBATCH --nodelist=n2
#SBATCH --partition=node2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --output=mergeMzml.sh.out

MZML_ROOT="/data/Storage_DAS03/home3/menteeE/minsea/Proteomics/mzml"

echo "merge mzML (Start)"

find "$MZML_ROOT" -type d -name individual | while read INDIV_DIR; do
    PROPRIETARY_DIR=$(dirname "$INDIV_DIR")
    SAMPLE_DIR=$(dirname "$PROPRIETARY_DIR")
    SAMPLE_NAME=$(basename "$SAMPLE_DIR")

    # 예: NCC_N111T112_Proteome_KU_20150109→ N111T112_Proteome
    CORE_NAME=$(echo "$SAMPLE_NAME" | cut -d'_' -f2,3)

    OUT_FILE="$PROPRIETARY_DIR/${CORE_NAME}_merged.mzML"

    echo "[$SAMPLE_NAME] merging files in $INDIV_DIR"
    FileMerger -in "$INDIV_DIR"/*.mzML -out "$OUT_FILE"
    echo "[$SAMPLE_NAME] merged to $OUT_FILE"
done

echo "mzML merge (Finish)"