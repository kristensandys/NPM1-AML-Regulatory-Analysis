# Select the top 16 KEGG pathways for all five NPM1 subtype projects.
# Process the upregulated and downregulated results separately.
# Rank pathways by adjusted p-value, followed by raw p-value.

top_n <- 16
top_label <- paste0("top", top_n)
kegg_library <- "KEGG_2026"

# 1. Define the five project locations
base_dir <- "/Users/kristensandys-thomas/Desktop/H3K27ac_diff_peak/h3k27ac_differential_peak_analysis"

project_dirs <- file.path(
  base_dir,
  c(
    "59_NPM1.1_differential_peak_project",
    "61_NPM1.2_differential_peak_project",
    "65_NPM1.3_differential_peak_project",
    "64_NPM1.4_differential_peak_project",
    "66_NPM1.5_differential_peak_project"
  )
)

# Label each project path with its NPM1 subtype.
names(project_dirs) <- paste0("NPM1.", 1:5)

# 2. Process the two KEGG enrichment tables from every NPM1 project
for (subtype in names(project_dirs)) {
  project_dir <- project_dirs[[subtype]]

  input_dir <- file.path(project_dir, "KEGG")
  output_dir <- file.path(
    input_dir,
    paste0(top_label, "_KEGG_pathways_3kb_group_specific")
  )

  dir.create(output_dir, showWarnings = FALSE)

  # define the input files 
  input_files <- file.path(
    input_dir,
    c(
      paste0(
        subtype,
        "_UPREGULATED_unique_3kb_GREAT_candidate_genes_",
        kegg_library,
        "_enrichment.tsv"
      ),
      paste0(
        subtype,
        "_DOWNREGULATED_unique_3kb_GREAT_candidate_genes_",
        kegg_library,
        "_enrichment.tsv"
      )
    )
  )

  for (input_file in input_files) {
    kegg <- read.delim(input_file, check.names = FALSE)

    # Treat the p-value columns as numbers.
    kegg$p_value <- as.numeric(kegg$p_value)
    kegg$adjusted_p_value <- as.numeric(kegg$adjusted_p_value)

    # Rank by adjusted p-value. If these are equal, rank by raw p-value.
    kegg <- kegg[
      order(kegg$adjusted_p_value, kegg$p_value, na.last = TRUE),
    ]

    # Keep the first 16 pathways after ranking.
    top_pathways <- head(kegg, top_n)

    output_name <- sub(
      "_enrichment\\.tsv$",
      paste0("_", top_label, ".tsv"),
      basename(input_file)
    )

    output_file <- file.path(output_dir, output_name)

    write.table(
      top_pathways,
      output_file,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      na = "NA"
    )
  }
}
