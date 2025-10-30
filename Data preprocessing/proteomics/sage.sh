#!/bin/bash
#SBATCH -J minsea_global_sage_start
#SBATCH --nodelist=n2
#SBATCH --partition=node2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=.global_sage.out

ORIG_JSON="./JSON/globalProteomics.json"
PATH_TXT="./JSON/mzmlPath.txt"
BATCH_SIZE=5
OUTPUT_DIR="./JSON"
OUTPUT_SAGE_DIR="./sage_res"

read -r line < "$PATH_TXT"

IFS=',' read -ra mzml_paths <<< "$line"

for i in "${!mzml_paths[@]}"; do
  mzml_paths[$i]="${mzml_paths[$i]//\"/}"
done

TOTAL=${#mzml_paths[@]}

count=1
for ((i=0; i <TOTAL; i+=BATCH_SIZE)); do
        part_paths=("${mzml_paths[@]:i:BATCH_SIZE}")
        jq_array=$(printf '%s\n' "${part_paths[@]}" | jq -R . | jq -s .)

        idx=$(printf "%02d" $count)
        out_json="${OUTPUT_DIR}/global_Proteomics_part_${idx}.json"

        jq --argjson paths "$jq_array" '.mzml_paths = $paths' "$ORIG_JSON" > "$out_json"

        echo "Created: $out_json"

        sbatch <<EOT
#!/bin/bash
#SBATCH -J minsea_global_sage_${idx}
#SBATCH --nodelist=n2
#SBATCH --partition=node2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --output=.global_sage_${idx}_%j.out

echo "sage proteomics part $idx (Start)"

sage $out_json -o ${OUTPUT_SAGE_DIR}/sage_res_${idx}

echo "sage proteomics part $idx (End)"
EOT

        ((count++))
done