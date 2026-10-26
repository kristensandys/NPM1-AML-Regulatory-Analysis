# Extract genes associated with peaks from the NPM1.2 versus NPM1.1 comparison.
# Process the upregulated and downregulated GREAT results separately.
# If the same gene is linked to multiple peaks, keep the peak located closest
# to the gene's transcription start site.

# 1. define the directory you will run the script from 
project_dir <- "/Users/kristensandys-thomas/Desktop/H3K27ac_diff_peak/test"

input_dir <- file.path(project_dir, "03_great", "raw")
output_dir <- file.path(project_dir, "03_great", "cleaned")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# define the paths to the upregulated and downregulated GREAT mapping files produced earlier.  
input_files <- file.path(
  input_dir,
  c(
    paste0(
      "NPM1.2_vs_NPM1.1_NPM1.2_UPREGULATED_regions_top100_",
      "GREAT_region_gene_mapping.tsv"
    ),
    paste0(
      "NPM1.2_vs_NPM1.1_NPM1.2_DOWNREGULATED_regions_top100_",
      "GREAT_region_gene_mapping.tsv"
    )
  )
)

# 2. Process the upregulated and downregulated GREAT mapping files
for (input_file in input_files) {
  
  mapping <- read.delim(
    input_file,
    check.names = FALSE
  )
  # Remove peaks without an associated gene.
  mapping <- mapping[!is.na(mapping$associated_genes) & mapping$associated_genes != "", ]

  # Separate semicolon-delimited genes and TSS distances into lists.
  gene_list <- strsplit(mapping$associated_genes,";", fixed = TRUE)

  distance_list <- strsplit(as.character(mapping$distance_to_tss_bp),";", fixed = TRUE)

  # Confirm that every associated gene has one corresponding TSS distance.
  if (any(lengths(gene_list) != lengths(distance_list))) {
    stop(
      "The number of genes and distances does not match in: ",
      input_file
    )
  }
  # Convert the lists into paired gene and numeric TSS distance vectors.
  genes <- trimws(unlist(gene_list, use.names = FALSE))
  distances <- trimws(unlist(distance_list, use.names = FALSE))
  distances <- as.numeric(distances)

  cleaned <- data.frame(
    gene = genes,
    distance_to_tss_bp = distances,
    absolute_distance_to_tss_bp = abs(distances)
  )
 
  # Sort by gene and absolute distance, the closest association appears first for every gene.
  cleaned <- cleaned[
    order(cleaned$gene, cleaned$absolute_distance_to_tss_bp), 
  ]

  # Keep each gene once, retaining its closest peak-to-TSS association.
  cleaned <- cleaned[!duplicated(cleaned$gene), ]

  # order the final table alphabetically by gene name 
  cleaned <- cleaned[order(cleaned$gene), ]

  rownames(cleaned) <- NULL

  # 3. Define and write the output files
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
