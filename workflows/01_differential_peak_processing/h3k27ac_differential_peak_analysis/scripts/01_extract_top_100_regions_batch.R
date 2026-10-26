# Stage 1: extract the top 100 differential peaks in both directions from every DESeq2 result file.

# Positive and negative log2FoldChange values are interpreted relative to the
# target NPM1 subtype named first in each comparison folder. 


top_n <- 100

# 1. Define the project locations
base_dir <- paste0(
  "/Users/kristensandys-thomas/Desktop/H3K27ac_diff_peak/",
  "h3k27ac_differential_peak_analysis"
)

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

# 2. Defined functions used for every comparison
# Write the genomic peak coordinates as row names in each output table.
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

# DESeq2 reports the first condition against the second. Return 1 when the
# target is first and -1 when the target is second, so every result can be
# described relative to the target subtype.
get_target_sign <- function(target_subtype, deseq_file) {
  
  #change, for example,NPM1.2 to NPM1S2 as this is the naming convention in the DESeq2 files. 
  target_code <- sub("\\.", "S", target_subtype)
  filename <- basename(deseq_file)
  
  # extract the two NPM1 groups from the filename. 
  contrast <- regmatches(
    filename,
    regexec(
      "condition_(NPM1S[1-5])_(NPM1S[1-5])\\.diffexp\\.tsv$",
      filename
    )
  )[[1]]
  
  # keep the origional fold-change when the target group is first 
  if (target_code == contrast[2]) return(1)
  
  # reverse the fold-change when the target is second. 
  if (target_code == contrast[3]) return(-1)
}

# 3. Process each NPM1 subtype project
for (project_dir in project_dirs) {
  
  deseq_files <- list.files(
    project_dir,
    pattern = "\\.diffexp\\.tsv$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  # keep only files located in a comparison's "01_deseq2" folder. 
  deseq_files <- sort(
    deseq_files[grepl("/comparisons/[^/]+/01_deseq2/", deseq_files)]
  )
  
  # keep and store counts 
  peak_counts <- data.frame()
  selected_up <- character()
  selected_down <- character()
  
  # 4. Process each DESeq2 file in the current project
  for (deseq_file in deseq_files) {
    
    comparison_dir <- dirname(dirname(deseq_file))
    comparison <- basename(comparison_dir)
    
    # First subtype in the comparison folder is the target subtype
    target_subtype <- sub("_vs_.*$", "", comparison)
    target_sign <- get_target_sign(target_subtype, deseq_file)
    
    peaks <- read.delim(
      deseq_file,
      row.names = 1,
      check.names = FALSE
    )
    
    peaks$region <- rownames(peaks)

    # remove peaks with "NA" values under log2FC
    peaks <- peaks[!is.na(peaks$log2FoldChange), ]
    
    # 5. Separate peaks into the two directions
    # Orient the fold change to the target subtype.
    target_log2fc <- peaks$log2FoldChange * target_sign
    
    up_index <- target_log2fc > 0
    down_index <- target_log2fc < 0
    
    up_peaks <- peaks[up_index, ]
    down_peaks <- peaks[down_index, ]
    
    # 6. Rank the peaks and select the top 100
    # Largest target-oriented positive fold changes first.
    up_peaks <- up_peaks[
      order(target_log2fc[up_index], decreasing = TRUE),
    ]
    
    # Most negative target-oriented fold changes first.
    down_peaks <- down_peaks[
      order(target_log2fc[down_index]),
    ]
    
    number_up_found <- nrow(up_peaks)
    number_down_found <- nrow(down_peaks)
    
    up_peaks <- head(up_peaks, top_n)
    down_peaks <- head(down_peaks, top_n)
    
    selected_up <- c(selected_up, up_peaks$region)
    selected_down <- c(selected_down, down_peaks$region)
    
    peak_counts <- rbind(
      peak_counts,
      data.frame(
        NPM_subtype = target_subtype,
        Comparison = comparison,
        Direction = c("UPREGULATED", "DOWNREGULATED"),
        Peaks_found = c(number_up_found, number_down_found),
        Peaks_selected = c(nrow(up_peaks), nrow(down_peaks))
      )
    )
    
    # 7. Write the two top-100 tables into the new directory
    output_dir <- file.path(comparison_dir, "02_top_regions")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    file_prefix <- paste0(comparison, "_", target_subtype)
    
    up_file <- file.path(
      output_dir,
      sprintf("%s_UPREGULATED_regions_top%s.tsv", file_prefix, top_n)
    )
    
    down_file <- file.path(
      output_dir,
      sprintf("%s_DOWNREGULATED_regions_top%s.tsv", file_prefix, top_n)
    )
    
    write_peak_table(up_peaks, up_file)
    write_peak_table(down_peaks, down_file)
  }
  
  # 8. Write summary tables for the current NPM1 subtype
  subtype <- unique(peak_counts$NPM_subtype)
  
  #create a directroy to stores these summaries 
  summary_dir <- file.path(project_dir, "summary")
  dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
  
  comparison_file <- file.path(
    summary_dir,
    paste0(subtype, "_peak_counts_by_comparison.tsv")
  )
  
  write.table(
    peak_counts,
    comparison_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  up_counts <- peak_counts[peak_counts$Direction == "UPREGULATED", ]
  down_counts <- peak_counts[peak_counts$Direction == "DOWNREGULATED", ]
  
  subtype_summary <- data.frame(
    NPM_subtype = subtype,
    Direction = c("UPREGULATED", "DOWNREGULATED"),
    Number_of_comparisons = c(nrow(up_counts), nrow(down_counts)),
    Total_peaks_found = c(
      sum(up_counts$Peaks_found),
      sum(down_counts$Peaks_found)
    ),
    Total_peaks_selected = c(
      sum(up_counts$Peaks_selected),
      sum(down_counts$Peaks_selected)
    ),
    Unique_selected_peaks = c(
      length(unique(selected_up)),
      length(unique(selected_down))
    )
  )
  
  summary_file <- file.path(
    summary_dir,
    paste0(subtype, "_total_peak_counts_by_direction.tsv")
  )
  
  write.table(
    subtype_summary,
    summary_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}
