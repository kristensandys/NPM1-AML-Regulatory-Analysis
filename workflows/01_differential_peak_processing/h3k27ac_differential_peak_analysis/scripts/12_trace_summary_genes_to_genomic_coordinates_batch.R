# Trace pooled candidate genes from all five NPM1 projects back to their
# genomic regions. Process upregulated and downregulated genes separately.

maximum_distance_bp <- 3000
directions <- c("UPREGULATED", "DOWNREGULATED")

# 1. Define the production project locations
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

# 2. Separate paired GREAT genes and distances into individual rows
expand_great_mapping <- function(mapping, input_file) {
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

  # Repeat each peak row once for every associated gene.
  source_rows <- rep(seq_len(nrow(mapping)), lengths(gene_lists))
  expanded <- mapping[source_rows, , drop = FALSE]

  expanded$gene <- unlist(gene_lists, use.names = FALSE)
  expanded$distance_to_tss_bp <- as.numeric(
    unlist(distance_lists, use.names = FALSE)
  )
  expanded$absolute_distance_to_tss_bp <- abs(expanded$distance_to_tss_bp)

  expanded$associated_genes <- NULL
  rownames(expanded) <- NULL
  expanded
}

# Store all subtype and direction trace tables.
all_traces <- list()

# 3. Process each NPM1 subtype project
for (subtype in names(project_dirs)) {
  project_dir <- project_dirs[[subtype]]
  summary_dir <- file.path(project_dir, "summary")

  # Trace the upregulated and downregulated genes separately.
  for (direction in directions) {
    prefix <- paste0(
      subtype, "_", direction, "_unique_3kb_GREAT_candidate_genes"
    )

    gene_file <- file.path(summary_dir, paste0(prefix, ".txt"))
    final_genes <- readLines(gene_file, warn = FALSE)

    raw_files <- list.files(
      file.path(project_dir, "comparisons"),
      pattern = paste0(
        "_", direction,
        "_regions_top[0-9]+_GREAT_region_gene_mapping\\.tsv$"
      ),
      recursive = TRUE,
      full.names = TRUE
    )

    raw_files <- sort(raw_files[
      grepl("/comparisons/[^/]+/03_great/raw/", raw_files)
    ])

    traced_parts <- list()

    # 4. Trace genes through each available comparison
    for (input_file in raw_files) {
      comparison_dir <- dirname(dirname(dirname(input_file)))
      comparison <- basename(comparison_dir)

      comparison_gene_file <- file.path(
        comparison_dir,
        "03_great",
        "cleaned_3kb",
        paste0(
          comparison, "_", subtype, "_", direction,
          "_GREAT_genes_within_3kb_only.txt"
        )
      )

      comparison_genes <- readLines(comparison_gene_file, warn = FALSE)
      trace_genes <- intersect(final_genes, comparison_genes)

      mapping <- read.delim(
        input_file, check.names = FALSE, stringsAsFactors = FALSE
      )

      expanded <- expand_great_mapping(mapping, input_file)

      # Retain selected genes whose peak is within 3 kb of the TSS.
      keep <- expanded$gene %in% trace_genes &
        !is.na(expanded$absolute_distance_to_tss_bp) &
        expanded$absolute_distance_to_tss_bp <= maximum_distance_bp

      traced <- expanded[keep, , drop = FALSE]

      traced$NPM_subtype <- subtype
      traced$Direction <- direction
      traced$Comparison <- comparison

      first_columns <- c(
        "NPM_subtype", "Direction", "Comparison", "gene", "peak_id",
        "chromosome", "start", "end", "distance_to_tss_bp",
        "absolute_distance_to_tss_bp"
      )

      traced <- traced[, c(
        first_columns, setdiff(names(traced), first_columns)
      ), drop = FALSE]

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

    all_traces[[length(all_traces) + 1L]] <- traced_all
  }
}

# 5. Combine all subtype and direction trace tables
combined_output <- file.path(
  base_dir,
  "ALL_NPM_unique_3kb_GREAT_candidate_genes_with_genomic_coordinates.tsv"
)

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
