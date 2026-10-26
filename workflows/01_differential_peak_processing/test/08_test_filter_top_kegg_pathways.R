# Select the top 16 KEGG pathways from the NPM1.2 enrichment results.
# Process the upregulated and downregulated results separately.
# Rank pathways by adjusted p-value, followed by raw p-value.

top_n <- 16
top_label <- paste0("top", top_n)
kegg_library <- "KEGG_2026"

# 1. input and output locations 
project_dir <- "/Users/kristensandys-thomas/Desktop/H3K27ac_diff_peak/test"

input_dir <- file.path(project_dir, "KEGG")
output_dir <- file.path(
  input_dir,
  paste0(top_label, "_KEGG_pathways_3kb_group_specific")
)

dir.create(output_dir, showWarnings = FALSE)

# define the two KEGG enrichment tables used as input. 
input_files <- file.path(
  input_dir,
  c(
    paste0(
      "NPM1.2_UPREGULATED_unique_3kb_GREAT_candidate_genes_",
      kegg_library,
      "_enrichment.tsv"
    ),
    paste0(
      "NPM1.2_DOWNREGULATED_unique_3kb_GREAT_candidate_genes_",
      kegg_library,
      "_enrichment.tsv"
    )
  )
)

# 2. Read each KEGG enrichment table into R
for (input_file in input_files) {
  kegg <- read.delim(input_file, check.names = FALSE)

  # Treat the p-value columns as numbers.
  kegg$p_value <- as.numeric(kegg$p_value)
  kegg$adjusted_p_value <- as.numeric(kegg$adjusted_p_value)

  #Rank by adjust p-value; if they are the same; then break ties with 
  #raw p-values
  kegg <- kegg[
    order(kegg$adjusted_p_value, kegg$p_value, na.last = TRUE),
  ]

  # Keep the first 16 pathways after ranking.
  top_pathways <- head(kegg, top_n)

  # Replace the enrichment suffix with the top-16 suffix.
  output_name <- sub(
    "_enrichment\\.tsv$",
    paste0("_", top_label, ".tsv"),
    basename(input_file)
  )

  output_file <- file.path(output_dir, output_name)

  # Write the selected pathways to a TSV file.
  write.table(
    top_pathways,
    output_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )
}
