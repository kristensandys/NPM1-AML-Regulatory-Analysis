#!/usr/bin/env bash

set -euo pipefail

new_root="/scratch/kristensandyst/h3k27ac_differential_peak_analysis"

declare -A project_for_subtype=(
  ["NPM1.1"]="59_NPM1.1_differential_peak_project"
  ["NPM1.2"]="61_NPM1.2_differential_peak_project"
  ["NPM1.3"]="65_NPM1.3_differential_peak_project"
  ["NPM1.4"]="64_NPM1.4_differential_peak_project"
  ["NPM1.5"]="66_NPM1.5_differential_peak_project"
)

subtypes=(
  "NPM1.1"
  "NPM1.2"
  "NPM1.3"
  "NPM1.4"
  "NPM1.5"
)

mkdir -p "$new_root/scripts"

for source_subtype in "${subtypes[@]}"; do
  project_name="${project_for_subtype[$source_subtype]}"
  project_dir="$new_root/$project_name"

  mkdir -p "$project_dir/comparisons" "$project_dir/summary"

  for target_subtype in "${subtypes[@]}"; do
    if [[ "$source_subtype" == "$target_subtype" ]]; then
      continue
    fi

    comparison="${source_subtype}_vs_${target_subtype}"
    comparison_dir="$project_dir/comparisons/$comparison"

    mkdir -p "$comparison_dir/01_deseq2"
  done
done

echo "Created the project structure under:"
echo "$new_root"
