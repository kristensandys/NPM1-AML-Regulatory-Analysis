
# Extract genes associated with peaks from allNPM1 comparisons.
# Process the upregulated and downregulated GREAT results separately.
# If the same gene is linked to multiple peaks, keep the peak located closest
# to the gene's transcription start site.

# define the project location 
base_dir <- "/Users/kristensandys-thomas/Desktop/H3K27ac_diff_peak/h3k27ac_differential_peak_analysis"

project_dirs <- file.path(
  base_dir,
  c(
    "59_NPM1.1_differential_peak_project",
    "61_NPM1.2_differential_peak_project",
    "64_NPM1.4_differential_peak_project",
    "65_NPM1.3_differential_peak_project",
    "66_NPM1.5_differential_peak_project"
  )
)

# find GREAT mapping files
input_files <- character()

for (project_dir in project_dirs) {
  project_files <- list.files(
    project_dir,
    pattern = "_GREAT_region_gene_mapping\\.tsv$",
    recursive = TRUE,
    full.names = TRUE
  )

  # Keep only files inside a comparison's 03_great/raw directory.
  project_files <- project_files[
    grepl("/comparisons/[^/]+/03_great/raw/", project_files)
  ]

  input_files <- c(input_files, project_files)
}

# 3. Process every available GREAT mapping file
for (input_file in input_files) {
  mapping <- read.delim(
    input_file,
    check.names = FALSE
  )

  # Remove peaks for which GREAT found no associated gene..
  mapping <- mapping[
    !is.na(mapping$associated_genes) &
      mapping$associated_genes != "",
  ]

  # Separate semicolon-delimited genes and TSS distances into lists.
  gene_list <- strsplit(
    mapping$associated_genes,
    ";",
    fixed = TRUE
  )

  distance_list <- strsplit(
    as.character(mapping$distance_to_tss_bp),
    ";",
    fixed = TRUE
  )

  # Confirm that every associated gene has one corresponding TSS distance.
  if (any(lengths(gene_list) != lengths(distance_list))) {
    stop(
      "The number of genes and distances does not match in: ",
      input_file
    )
  }

  # Convert the lists into paired gene and numeric TSS-distance vectors.
  genes <- trimws(unlist(gene_list, use.names = FALSE))
  distances <- trimws(unlist(distance_list, use.names = FALSE))
  distances <- as.numeric(distances)

  cleaned <- data.frame(
    gene = genes,
    distance_to_tss_bp = distances,
    absolute_distance_to_tss_bp = abs(distances)
  )

  # Sort by gene and absolute distance so that the closest association appears
  # first for every gene.
  cleaned <- cleaned[
    order(cleaned$gene, cleaned$absolute_distance_to_tss_bp),
  ]

  # Keep each gene once, retaining its closest peak-to-TSS association.
  cleaned <- cleaned[!duplicated(cleaned$gene), ]

  # order tables alphabetically by gene name 
  cleaned <- cleaned[order(cleaned$gene), ]
  rownames(cleaned) <- NULL

  # 4. Define the output location for the current comparison
  comparison_dir <- dirname(dirname(dirname(input_file)))
  output_dir <- file.path(comparison_dir, "03_great", "cleaned")

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # naming convention 
  prefix <- sub(
    "_regions_top[0-9]+_GREAT_region_gene_mapping\\.tsv$",
    "_GREAT_genes",
    basename(input_file)
  )

  distance_file <- file.path(
    output_dir,
    paste0(prefix, "_with_distance.tsv")
  )

  genes_file <- file.path(
    output_dir,
    paste0(prefix, "_only.txt")
  )

  # 5. Write the distance table and gene-only list
  write.table(
    cleaned,
    distance_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )

  writeLines(cleaned$gene, genes_file)
}
