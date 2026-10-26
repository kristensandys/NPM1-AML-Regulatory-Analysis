
# Trace the pooled NPM1.2 candidate genes back to their genomic regions.
# Process the upregulated and downregulated gene lists separately.


maximum_distance_bp <- 3000
subtype <- "NPM1.2"
directions <- c("UPREGULATED", "DOWNREGULATED")

# 1. Define the test input and output locations
project_dir <- "/Users/kristensandys-thomas/Desktop/H3K27ac_diff_peak/test"

summary_dir <- file.path(project_dir, "summary")
raw_dir <- file.path(project_dir, "03_great", "raw")
cleaned_3kb_dir <- file.path(project_dir, "03_great", "cleaned_3kb")

# 4. Separate paired GREAT genes and distances into individual rows
expand_great_mapping <- function(mapping, input_file) {
  # Remove peaks without an associated gene.
  mapping <- mapping[
    !is.na(mapping$associated_genes) & mapping$associated_genes != "",
    ,
    drop = FALSE
  ]

  gene_lists <- strsplit(as.character(mapping$associated_genes), ";", fixed = TRUE)

  distance_lists <- strsplit(
    as.character(mapping$distance_to_tss_bp), ";", fixed = TRUE
  )

  # Confirm that every gene remains paired with one TSS distance.
  if (any(lengths(gene_lists) != lengths(distance_lists))) {
    stop("Gene/distance pairing is inconsistent in: ", input_file)
  }

  # repeat each peak row once for every gene associated with that peak; and then
  # create one copy of the peak information for each associated gene. 
  source_rows <- rep(seq_len(nrow(mapping)), lengths(gene_lists))
  expanded <- mapping[source_rows, , drop = FALSE]

  # Place one associated gene in each row.
  expanded$gene <- unlist(gene_lists, use.names = FALSE)

  # Place the corresponding TSS distance in the same row as its gene.
  expanded$distance_to_tss_bp <- as.numeric(
    unlist(distance_lists, use.names = FALSE)
  )
  
  #calculate the distance 
  expanded$absolute_distance_to_tss_bp <- abs(expanded$distance_to_tss_bp)

  # Replace the packed gene column with one gene per row.
  expanded$associated_genes <- NULL
  rownames(expanded) <- NULL

  expanded
}

# create an empty list to store the co-ordiante trace 
all_traces <- list()

# 5. Trace the upregulated and downregulated genes separately
for (direction in directions) {
  prefix <- paste0(
    subtype, "_", direction, "_unique_3kb_GREAT_candidate_genes"
  )

  gene_file <- file.path(summary_dir, paste0(prefix, ".txt"))
  final_genes <- readLines(gene_file, warn = FALSE)

  raw_files <- list.files(
    raw_dir,
    pattern = paste0(
      "_", direction,
      "_regions_top[0-9]+_GREAT_region_gene_mapping\\.tsv$"
    ),
    full.names = TRUE
  )

  raw_files <- sort(raw_files)
  traced_parts <- list()

  for (input_file in raw_files) {
    filename <- basename(input_file)

    # Extract the comparison name from the raw GREAT filename.
    comparison <- strsplit(
      filename, paste0("_", subtype, "_"), fixed = TRUE
    )[[1]][1]

    comparison_gene_file <- file.path(cleaned_3kb_dir, paste0(
      comparison, "_", subtype, "_", direction,
      "_GREAT_genes_within_3kb_only.txt"
    ))

    # identify which pooled genes came from the current comparison 
    comparison_genes <- readLines(comparison_gene_file, warn = FALSE)
    # keep pooled genes that were found in this comaprison. 
    trace_genes <- intersect(final_genes, comparison_genes)
    
    #read the raw GREAT mapping .tsv into R
    mapping <- read.delim(input_file, check.names = FALSE, stringsAsFactors = FALSE)

    expanded <- expand_great_mapping(mapping, input_file)
   
    # Retain selected genes whose peak is within 3 kb of the TSS.
    keep <- expanded$gene %in% trace_genes &
      !is.na(expanded$absolute_distance_to_tss_bp) &
      expanded$absolute_distance_to_tss_bp <= maximum_distance_bp

    traced <- expanded[keep, , drop = FALSE]

    traced$NPM_subtype <- subtype
    traced$Direction <- direction
    traced$Comparison <- comparison
    
    # arrange the columns 
    first_columns <- c(
      "NPM_subtype", "Direction", "Comparison", "gene", "peak_id",
      "chromosome", "start", "end", "distance_to_tss_bp",
      "absolute_distance_to_tss_bp"
    )
    
    # define which columns appear first 
    traced <- traced[, c(
      first_columns, setdiff(names(traced), first_columns)
    ), drop = FALSE]
    
    # store each comparison 
    traced_parts[[length(traced_parts) + 1L]] <- traced
  }

  
    traced_all <- do.call(rbind, traced_parts)
    traced_all <- unique(traced_all)
    traced_all <- traced_all[
      order(traced_all$gene, traced_all$Comparison,
        traced_all$absolute_distance_to_tss_bp),
    ]
    rownames(traced_all) <- NULL

  output_file <- file.path(
    summary_dir, paste0(prefix, "_with_genomic_coordinates.tsv")
  )

  write.table(
    traced_all,
    output_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )

  # Store this direction for the final combined table.
  all_traces[[length(all_traces) + 1L]] <- traced_all
}

 
# 6. Combine the upregulated and downregulated trace tables
combined_output <- file.path(
  project_dir,
  "TEST_NPM1.2_unique_3kb_GREAT_candidate_genes_with_genomic_coordinates.tsv"
)

# Combine the upregulated and downregulated trace tables.
combined <- do.call(rbind, all_traces)
rownames(combined) <- NULL

write.table(
  combined,
  combined_output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)
