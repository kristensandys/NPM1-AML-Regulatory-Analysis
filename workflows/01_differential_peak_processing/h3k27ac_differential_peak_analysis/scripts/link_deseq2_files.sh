#!/usr/bin/env bash

# Link each pairwise DESeq2 result into both directed NPM1 comparisons.
# The source files remain unchanged; only symbolic links are created.

set -euo pipefail

new_root="/scratch/kristensandyst/h3k27ac_differential_peak_analysis"
old_root="/scratch/kristensandyst/diff_peak_analysis"

declare -A project_for_subtype=(
  ["NPM1.1"]="59_NPM1.1_differential_peak_project"
  ["NPM1.2"]="61_NPM1.2_differential_peak_project"
  ["NPM1.3"]="65_NPM1.3_differential_peak_project"
  ["NPM1.4"]="64_NPM1.4_differential_peak_project"
  ["NPM1.5"]="66_NPM1.5_differential_peak_project"
)

sources=(
  "$old_root/47_NPM1.3_NPM1.1/results/deseq2/GRCh38-condition_NPM1S3_NPM1S1.diffexp.tsv"
  "$old_root/48_NPM1.5_NPM1.1/results/deseq2/GRCh38-condition_NPM1S5_NPM1S1.diffexp.tsv"
  "$old_root/49_NPM1.2_NPM1.3/results/deseq2/GRCh38-condition_NPM1S3_NPM1S2.diffexp.tsv"
  "$old_root/50_NPM1.4_NPM1.2/results/deseq2/GRCh38-condition_NPM1S4_NPM1S2.diffexp.tsv"
  "$old_root/51_NPM1.5_NPM1.2/results/deseq2/GRCh38-condition_NPM1S5_NPM1S2.diffexp.tsv"
  "$old_root/52_NPM1.4_NPM1.3/results/deseq2/GRCh38-condition_NPM1S4_NPM1S3.diffexp.tsv"
  "$old_root/53_NPM1.5_NPM1.3/results/deseq2/GRCh38-condition_NPM1S5_NPM1S3.diffexp.tsv"
  "$old_root/55_NPM1.2_NPM1.1/results/deseq2/GRCh38-condition_NPM1S2_NPM1S1.diffexp.tsv"
  "$old_root/56_NPM1.1_NPM1.4/results/deseq2/GRCh38-condition_NPM1S4_NPM1S1.diffexp.tsv"
)

linked=0
skipped=0
errors=0

for source_file in "${sources[@]}"; do
  if [[ ! -f "$source_file" ]]; then
    echo "MISSING SOURCE: $source_file" >&2
    ((errors += 1))
    continue
  fi

  filename="$(basename "$source_file")"

  if [[ "$filename" =~ NPM1S([1-5])_NPM1S([1-5]) ]]; then
    first="NPM1.${BASH_REMATCH[1]}"
    second="NPM1.${BASH_REMATCH[2]}"
  else
    echo "UNRECOGNIZED FILENAME: $source_file" >&2
    ((errors += 1))
    continue
  fi

  comparisons=(
    "${first}_vs_${second}"
    "${second}_vs_${first}"
  )

  for comparison in "${comparisons[@]}"; do
    project_subtype="${comparison%%_vs_*}"
    project_name="${project_for_subtype[$project_subtype]}"
    destination_dir="$new_root/$project_name/comparisons/$comparison/01_deseq2"
    destination="$destination_dir/$filename"

    if [[ ! -d "$destination_dir" ]]; then
      echo "MISSING DESTINATION: $destination_dir" >&2
      ((errors += 1))
      continue
    fi

    if [[ -e "$destination" || -L "$destination" ]]; then
      echo "SKIP EXISTING: $destination"
      ((skipped += 1))
      continue
    fi

    ln -s -- "$source_file" "$destination"
    echo "LINKED: $destination -> $source_file"
    ((linked += 1))
  done
done

echo
echo "Linked:  $linked"
echo "Skipped: $skipped"
echo "Errors:  $errors"
echo

echo "DESeq2 links by destination:"
find "$new_root" \
  -type l \
  -path '*/comparisons/*/01_deseq2/*.diffexp.tsv' \
  -exec ls -l -- {} \;

link_count="$({
  find "$new_root" \
    -type l \
    -path '*/comparisons/*/01_deseq2/*.diffexp.tsv'
} | wc -l)"

echo
echo "Total DESeq2 links: $link_count"

echo
echo "Comparison folders without a DESeq2 link:"
empty_count=0

while IFS= read -r destination_dir; do
  count="$({
    find "$destination_dir" -maxdepth 1 -type l -name '*.diffexp.tsv'
  } | wc -l)"

  if [[ "$count" -eq 0 ]]; then
    echo "$destination_dir"
    ((empty_count += 1))
  fi
done < <(
  find "$new_root" -type d -path '*/comparisons/*/01_deseq2' | sort
)

if [[ "$empty_count" -eq 0 ]]; then
  echo "None"
fi

echo
echo "Broken links:"
broken_count=0

while IFS= read -r link; do
  if [[ ! -e "$link" ]]; then
    echo "$link"
    ((broken_count += 1))
  fi
done < <(
  find "$new_root" \
    -type l \
    -path '*/comparisons/*/01_deseq2/*.diffexp.tsv' |
    sort
)

if [[ "$broken_count" -eq 0 ]]; then
  echo "None"
fi

if [[ "$errors" -gt 0 || "$broken_count" -gt 0 ]]; then
  exit 1
fi
