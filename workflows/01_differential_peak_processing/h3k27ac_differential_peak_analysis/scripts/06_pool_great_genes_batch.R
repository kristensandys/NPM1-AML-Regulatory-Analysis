# Pool the GREAT-associated genes located within 3 kb of a TSS for all five
# NPM1 subtype projects.
# Process the upregulated and downregulated gene lists separately.

# 1. project location
base_dir <- "/Users/kristensandys-thomas/Desktop/H3K27ac_diff_peak/h3k27ac_differential_peak_analysis"


project_dirs <- file.path(
  base_dir,
  c(
    "NPM1.1" = "59_NPM1.1_differential_peak_project",
    "NPM1.2" = "61_NPM1.2_differential_peak_project",
    "NPM1.3" = "65_NPM1.3_differential_peak_project",
    "NPM1.4" = "64_NPM1.4_differential_peak_project",
    "NPM1.5" = "66_NPM1.5_differential_peak_project"
  )
)

# Label each project path with its NPM1 subtype for the loop below.
names(project_dirs) <- paste0("NPM1.", 1:5)

# 2. Pool the gene lists for one subtype and direction
pool_direction <- function(project_dir, subtype, direction) {
  
  # search for gene lists inside the project's comparison folders
  input_dir <- file.path(project_dir, "comparisons")
  output_dir <- file.path(project_dir, "summary")


  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  great_pattern <- paste0(
    "_", direction,
    "_GREAT_genes_within_3kb_only\\.txt$"
  )
  
# find the input files (the GREAT files)
  great_files <- list.files(
    input_dir,
    pattern = great_pattern,
    recursive = TRUE,
    full.names = TRUE
  )

  # Identify the comparison represented by each gene list.
  comparison_names <- basename(
    dirname(dirname(dirname(great_files)))
  )

  # Read each comparison's gene list.
  gene_lists <- lapply(
    great_files,
    readLines,
    warn = FALSE
  )

  names(gene_lists) <- comparison_names

  # Combine all genes, remove duplicates, and order them alphabetically.
  pooled_genes <- sort(unique(unlist(gene_lists, use.names = FALSE)))

  # 3. Define the output filenames
  prefix <- paste0(
    subtype,
    "_",
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

  # Write the pooled unique gene list.
  writeLines(pooled_genes, pooled_file)

  # 4. Record whether each gene is present in each comparison
  presence_table <- data.frame(Gene = pooled_genes)

  for (comparison in comparison_names) {
    presence_table[[comparison]] <- ifelse(
      pooled_genes %in% gene_lists[[comparison]],
      "YES",
      "NO"
    )
  }

  # Count the number of comparisons containing each gene.
  presence_table$Frequency <- rowSums(
    presence_table[comparison_names] == "YES"
  )

  # Show the most frequently observed genes first, then order ties by gene.
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

# 5. Pool both directions for every NPM1 subtype
for (subtype in names(project_dirs)) {
  project_dir <- project_dirs[[subtype]]
  pool_direction(project_dir, subtype, "UPREGULATED")
  pool_direction(project_dir, subtype, "DOWNREGULATED")
}
