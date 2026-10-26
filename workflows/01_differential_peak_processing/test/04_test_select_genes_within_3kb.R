# Select GREAT-associated genes located within 3 kb of TSS for the NPM1.2 versus NPM1.1 comparison.
# Process the upregulated and downregulated gene-distance tables separately.

maximum_distance_bp <- 3000

# 1. Define the input and output locations
project_dir <- "/Users/kristensandys-thomas/Desktop/H3K27ac_diff_peak/test"

input_dir <- file.path(project_dir, "03_great", "cleaned")
output_dir <- file.path(project_dir, "03_great", "cleaned_3kb")

# create the output directory 
dir.create(output_dir, showWarnings = FALSE)

# define the input files
input_files <- file.path(
  input_dir,
  c(
    paste0(
      "NPM1.2_vs_NPM1.1_NPM1.2_UPREGULATED_",
      "GREAT_genes_with_distance.tsv"
    ),
    paste0(
      "NPM1.2_vs_NPM1.1_NPM1.2_DOWNREGULATED_",
      "GREAT_genes_with_distance.tsv"
    )
  )
)

# 2. Process the upregulated and downregulated gene-distance tables
for (input_file in input_files) {
  genes <- read.delim(
    input_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  # Keep genes located no more than 3,000 bp from a TSS.
  keep <- (
    !is.na(genes$absolute_distance_to_tss_bp) &
      genes$absolute_distance_to_tss_bp <= maximum_distance_bp
  )

  selected <- genes[keep, , drop = FALSE]

  # Order the selected genes alphabetically for easier inspection.
  selected <- selected[order(selected$gene), , drop = FALSE]
  rownames(selected) <- NULL

  # 3. Define the Stage 4 output filenames
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

  # 4. Write the distance table and gene-only list
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
