
# Keep top-16 KEGG pathways with at least two overlapping genes.
# Identify pathways shared by at least two of the five NPM1 subtypes.

top_n <- 16
top_label <- paste0("top", top_n)
min_genes <- 2
min_groups <- 2
groups <- paste0("NPM1.", 1:5)

# 1. Define the production input directory
base_dir <- "/Users/kristensandys-thomas/Desktop/H3K27ac_diff_peak/h3k27ac_differential_peak_analysis"

input_dir <- file.path(
  base_dir,
  paste0("compared_", top_label, "_KEGG_pathways_3kb_group_specific")
)

input_files <- file.path(
  input_dir,
  c(
    paste0(
      "ALL_NPM_UPREGULATED_3kb_GREAT_KEGG_",
      top_label,
      "_combined_with_genes.tsv"
    ),
    paste0(
      "ALL_NPM_DOWNREGULATED_3kb_GREAT_KEGG_",
      top_label,
      "_combined_with_genes.tsv"
    )
  )
)

# 2. Count the unique genes overlapping one KEGG pathway
count_genes <- function(x) {
  genes <- trimws(strsplit(as.character(x), ";", fixed = TRUE)[[1]])
  length(unique(genes))
}

# 3. Process the upregulated and downregulated combined tables
for (input_file in input_files) {
  file_stem <- sub("_combined_with_genes\\.tsv$", "", input_file)

  backup_file <- paste0(
    file_stem,
    "_combined_with_genes_unfiltered.tsv"
  )

  shared_file <- paste0(file_stem, "_shared_summary.tsv")

  # Preserve the original combined table before applying the gene cutoff.
  if (!file.exists(backup_file)) {
    file.copy(input_file, backup_file)
  }

  combined <- read.delim(backup_file, check.names = FALSE)

  # Count the overlapping genes for every pathway row.
  combined$overlapping_gene_count <- vapply(
    combined$overlapping_genes,
    count_genes,
    integer(1)
  )

  # Keep pathways containing at least two overlapping genes.
  retained <- combined[
    combined$overlapping_gene_count >= min_genes,
  ]

  write.table(
    retained,
    input_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )

  # Keep one pathway-subtype pair before making the presence table.
  pathway_groups <- unique(retained[, c("term", "NPM_group")])
  all_pathways <- sort(unique(pathway_groups$term))
  presence_table <- data.frame(Pathway = all_pathways)

  # Record which NPM1 subtypes contain each pathway.
  for (subtype in groups) {
    subtype_pathways <- pathway_groups$term[
      pathway_groups$NPM_group == subtype
    ]

    presence_table[[subtype]] <- ifelse(
      presence_table$Pathway %in% subtype_pathways,
      "YES",
      "NO"
    )
  }

  presence_table$Frequency <- rowSums(
    presence_table[groups] == "YES"
  )

  presence_table <- presence_table[
    order(-presence_table$Frequency, presence_table$Pathway),
  ]

  rownames(presence_table) <- NULL

  # Keep pathways shared by at least two NPM1 subtypes.
  shared_table <- presence_table[
    presence_table$Frequency >= min_groups,
  ]

  write.table(
    shared_table,
    shared_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )
}
