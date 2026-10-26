# Combine the top-16 KEGG pathways across all five NPM1 subtypes.
# Process upregulated and downregulated results separately.

top_n <- 16
top_label <- paste0("top", top_n)
kegg_library <- "KEGG_2026"

groups <- paste0("NPM1.", 1:5)
directions <- c("UPREGULATED", "DOWNREGULATED")

# 1. project location
base_dir <- "/Users/kristensandys-thomas/Desktop/H3K27ac_diff_peak/h3k27ac_differential_peak_analysis"


output_dir <- file.path(
  base_dir,
  paste0("compared_", top_label, "_KEGG_pathways_3kb_group_specific")
)

# creates new directory to store resutls. 
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 2. Find the top-16 KEGG tables in the five subtype projects
top_suffix <- paste0("_", kegg_library, "_", top_label, "\\.tsv$")
top_folder <- paste0(
  "/KEGG/", top_label, "_KEGG_pathways_3kb_group_specific/"
)

all_files <- list.files(
  base_dir,
  pattern = top_suffix,
  recursive = TRUE,
  full.names = TRUE
)

# keep only files from the top-16 KEGG results folders. 
all_files <- all_files[
  grepl(top_folder, all_files, fixed = TRUE)
]

all_files <- sort(all_files)

# 3. Combine the subtype tables for each direction
for (direction in directions) {
  direction_pattern <- paste0(
    "_", direction,
    "_unique_3kb_GREAT_candidate_genes_",
    kegg_library,
    "_", top_label,
    "\\.tsv$"
  )

  direction_files <- all_files[
    grepl(direction_pattern, basename(all_files))
  ]
  
  # create an empty list to store each subtype's KEGG table. 
  subtype_tables <- list()

  for (input_file in direction_files) {
    filename <- basename(input_file)
    subtype <- sub("_.*$", "", filename)

    kegg <- read.delim(input_file, check.names = FALSE)
    # add columns to results 
    kegg$NPM_group <- subtype
    kegg$Direction <- direction

    subtype_tables[[subtype]] <- kegg
  }

  # Keep the NPM1 subtype order consistent.
  present_groups <- groups[groups %in% names(subtype_tables)]
  subtype_tables <- subtype_tables[present_groups]

  combined <- do.call(rbind, subtype_tables)
  rownames(combined) <- NULL

  output_name <- paste0(
    "ALL_NPM_", direction,
    "_3kb_GREAT_KEGG_", top_label,
    "_combined_with_genes.tsv"
  )

  output_file <- file.path(output_dir, output_name)

  write.table(
    combined,
    output_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )
}
