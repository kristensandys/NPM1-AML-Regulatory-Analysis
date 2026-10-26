# Pool the GREAT-associated genes located within 3 kb of a TSS for the
# NPM1.2 versus NPM1.1; NPM1.2 vs NPM1.3, and NPM1.2 s NPM1.4 comparisons.
# Process the upregulated and downregulated gene lists separately.

# 1. project location
project_dir <- "/Users/kristensandys-thomas/Desktop/H3K27ac_diff_peak/test"

input_dir <- file.path(project_dir, "03_great", "cleaned_3kb")
output_dir <- file.path(project_dir, "summary")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# 2. Pool the gene list for one direction
pool_direction <- function(direction) {
  great_pattern <- paste0(
    "_", direction,
    "_GREAT_genes_within_3kb_only\\.txt$"
  )

  great_files <- list.files(
    input_dir,
    pattern = great_pattern,
    full.names = TRUE
  )

  comparison_names <- sub(
    "_NPM1\\.2_.*$",
    "",
    basename(great_files)
  )

  # Read each comparison's gene list.
  gene_lists <- lapply(
    great_files,
    readLines,
    warn = FALSE
  )
  names(gene_lists) <- comparison_names

  pooled_genes <- sort(unique(unlist(gene_lists, use.names = FALSE)))

  prefix <- paste0(
    "NPM1.2_",
    direction,
    "_unique_3kb_GREAT_candidate_genes"
  )

  pooled_file <- file.path(
    output_dir,
    paste0(prefix, ".txt")
  )

  presence_file <- file.path(
    output_dir,
    paste0(prefix, "_presence_absence.tsv")
  )

  # Write the unique pooled gene list in alphabetical order.
  writeLines(pooled_genes, pooled_file)

  # Record whether each gene is present in each test comparison.
  presence_table <- data.frame(Gene = pooled_genes)

  for (comparison in comparison_names) {
    presence_table[[comparison]] <- ifelse(
      pooled_genes %in% gene_lists[[comparison]],
      "YES",
      "NO"
    )
  }

  presence_table$Frequency <- rowSums(
    presence_table[comparison_names] == "YES"
  )

  presence_table <- presence_table[
    order(-presence_table$Frequency, presence_table$Gene),
  ]

  rownames(presence_table) <- NULL

  write.table(
    presence_table,
    presence_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

# 4. Pool the upregulated and downregulated directions
pool_direction("UPREGULATED")
pool_direction("DOWNREGULATED")
