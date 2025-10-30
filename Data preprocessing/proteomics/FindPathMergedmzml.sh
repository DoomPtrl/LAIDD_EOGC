#!/usr/bin/env bash
  
MZML_ROOT="./mzml"
OUTPUT_TSV="MergedMzmlPath.tsv"

> "$OUTPUT_TSV"
declare -A grouped_paths

while IFS= read -r FILE; do
        TYPE=$(basename "$(dirname "$(dirname "$(dirname "$FILE")")")")
        FILE_QUOTED="\"$FILE\""

        if [ -z "${grouped_paths[$TYPE]}" ]; then
                grouped_paths[$TYPE]="$FILE_QUOTED"
        else
                grouped_paths[$TYPE]="${grouped_paths[$TYPE]},$FILE_QUOTED"
        fi
done < <(find "$MZML_ROOT" -type f -name "*_merged.mzML")

echo -e "type\tpaths" > "$OUTPUT_TSV"
for TYPE in "${!grouped_paths[@]}"; do
        echo -e "${TYPE}\t${grouped_paths[$TYPE]}" >> "$OUTPUT_TSV"
done

echo "Merged mzml tsv created path: $OUTPUT_TSV"