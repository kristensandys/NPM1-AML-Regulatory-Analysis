# Select GREAT-associated genes located within 3 kb of TSS from all NPM1 comparisons.
# Process the upregulated and downregulated gene-distance tables separately.

maximum_distance_bp <- 3000

# 1. project location
base_dir <- "/Users/kristensandys-thomas/Desktop/H3K27ac_diff_peak/h3k27ac_differential_peak_analysis"


# 2. find all gene-distance tables
input_files <- list.files(
  base_dir,
  pattern = "_GREAT_genes_with_distance\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

# 3. process gene-distance tables
for (input_file in input_files) {
  genes <- read.delim(
    input_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  # Keep genes located within 3,000 bp from a TSS.
  keep <- (
    !is.na(genes$absolute_distance_to_tss_bp) &
      genes$absolute_distance_to_tss_bp <= maximum_distance_bp
  )

  selected <- genes[keep, , drop = FALSE]

  # order genes alphabetically
  selected <- selected[order(selected$gene), , drop = FALSE]
  rownames(selected) <- NULL

  # 4. Define the output location for the current comparison
  comparison_dir <- dirname(dirname(dirname(input_file)))
  output_dir <- file.path(comparison_dir, "03_great", "cleaned_3kb")

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  prefix <- sub(
    "_with_distance\\.tsv$",
    "",
    basename(input_file)
  )

  distance_file <- file.path(
    output_dir,
    paste0(prefix, "_within_3kb_with_distance.tsv")
  )

  gene_file <- file.path(
    output_dir,
    paste0(prefix, "_within_3kb_only.txt")
  )

  # 5. Write the distance table and gene-only list
  write.table(
    selected,
    distance_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )

  writeLines(selected$gene, gene_file)
}
