#!/usr/bin/env Rscript

# Test script 01: extract the top 100 differential peaks in both directions
# from one DESeq2 result file.
#
# Test comparison: NPM1.2 versus NPM1.1
# Positive log2FoldChange = upregulated in NPM1.2
# Negative log2FoldChange = downregulated in NPM1.2
# Peaks are ranked by log2 fold change; padj is not used for selection.

top_n <- 100

# 1. Define the input and output locations
project_dir <- "/Users/kristensandys-thomas/Desktop/H3K27ac_diff_peak/test"

deseq_file <- file.path(
  project_dir,
  "deseq",
  "GRCh38-condition_NPM1S2_NPM1S1.diffexp.tsv"
)

output_dir <- file.path(project_dir, "02_top_regions")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


peaks <- read.delim(
  deseq_file,
  row.names = 1,
  check.names = FALSE
)

# genomic coordinates stored as row names.
peaks$region <- rownames(peaks)

# Remove peaks without a log2-fold-change value as some of my files have "NA" in this column near the end
peaks <- peaks[!is.na(peaks$log2FoldChange), ]

# Separate peaks into the two directions
up_peaks <- peaks[peaks$log2FoldChange > 0, ]
down_peaks <- peaks[peaks$log2FoldChange < 0, ]

number_up_found <- nrow(up_peaks)
number_down_found <- nrow(down_peaks)

# Rank the peaks and select the top 100
# Largest positive fold changes first.
up_peaks <- up_peaks[
  order(up_peaks$log2FoldChange, decreasing = TRUE),
]

# Most negative fold changes first.
down_peaks <- down_peaks[
  order(down_peaks$log2FoldChange),
]

up_peaks <- head(up_peaks, top_n)
down_peaks <- head(down_peaks, top_n)

# Write the two top-100 tables
write_peak_table <- function(peaks, output_file) {
  rownames(peaks) <- peaks$region
  peaks$region <- NULL

  write.table(
    peaks,
    output_file,
    sep = "\t",
    quote = FALSE,
    row.names = TRUE,
    col.names = NA,
    na = "NA"
  )
}

up_file <- file.path(
  output_dir,
  "NPM1.2_vs_NPM1.1_NPM1.2_UPREGULATED_regions_top100.tsv"
)

down_file <- file.path(
  output_dir,
  "NPM1.2_vs_NPM1.1_NPM1.2_DOWNREGULATED_regions_top100.tsv"
)

write_peak_table(up_peaks, up_file)
write_peak_table(down_peaks, down_file)


# make a summary table
peak_counts <- data.frame(
  NPM_subtype = "NPM1.2",
  Comparison = "NPM1.2_vs_NPM1.1",
  Direction = c("UPREGULATED", "DOWNREGULATED"),
  Peaks_found = c(number_up_found, number_down_found),
  Peaks_selected = c(nrow(up_peaks), nrow(down_peaks))
)

summary_file <- file.path(output_dir, "01_test_peak_counts.tsv")

write.table(
  peak_counts,
  summary_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
