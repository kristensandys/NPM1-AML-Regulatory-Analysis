# Combine top-16 KEGG pathway tables from the KEGG enrichment results. 
# Process the upregulated and downregulated results separately.

top_n <- 16
top_label <- paste0("top", top_n)
kegg_library <- "KEGG_2026"

groups <- paste0("NPM1.", 1:5)
directions <- c("UPREGULATED", "DOWNREGULATED")

# 1. the testing directory 
project_dir <- "/Users/kristensandys-thomas/Desktop/H3K27ac_diff_peak/test"

# define where the input files are: 
input_dir <- file.path(
  project_dir,
  "KEGG",
  paste0(top_label, "_KEGG_pathways_3kb_group_specific")
)

output_dir <- file.path(
  project_dir,
  paste0("compared_", top_label, "_KEGG_pathways_3kb_group_specific")
)

dir.create(output_dir, showWarnings = FALSE)

# 2. Find the top-16 KEGG tables in the test input directory
input_pattern <- paste0(
  "_", kegg_library,
  "_", top_label,
  "\\.tsv$"
)

all_files <- list.files(
  input_dir,
  pattern = input_pattern,
  full.names = TRUE
)

all_files <- sort(all_files)

# 3. Combine the available subtype tables for each direction
for (direction in directions) {
  direction_pattern <- paste0(
    "_", direction,
    "_unique_3kb_GREAT_candidate_genes_",
    kegg_library,
    "_", top_label,
    "\\.tsv$"
  )

  # keep files for the current direction 
  direction_files <- all_files[
    grepl(direction_pattern, basename(all_files))
  ]
  # create an empty list for sorting each subtype's KEGG table
  subtype_tables <- list()

  # Read each available subtype table and label its rows.
  for (input_file in direction_files) {
    filename <- basename(input_file)
    subtype <- sub("_.*$", "", filename)

    kegg <- read.delim(input_file, check.names = FALSE)
    
    # label results with its NPM1 subtype and direction. 
    kegg$NPM_group <- subtype
    kegg$Direction <- direction
    
    # store the table 
    subtype_tables[[subtype]] <- kegg
  }

  # Keep the NPM1 subtype order consistent.
  present_groups <- groups[groups %in% names(subtype_tables)]
  subtype_tables <- subtype_tables[present_groups]

  # Combine the subtype tables into one table.
  combined <- do.call(rbind, subtype_tables)
  rownames(combined) <- NULL

  output_name <- paste0(
    "TEST_ALL_NPM_",
    direction,
    "_3kb_GREAT_KEGG_",
    top_label,
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
