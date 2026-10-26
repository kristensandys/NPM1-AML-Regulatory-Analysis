# assign genes to selected differential peaks from NPM1.2 vs NPM1.1 comparison. 
# Both the upregulated and downregulated top-100 peak files are processed.
# Gene assignment and TSS distances are based on peak midpoints.
# Negative distances are upstream of the TSS and positive distances are
# downstream of the TSS.

library(rGREAT)
library(GenomicRanges)
library(IRanges)

# 1. Define the input and output locations
project_dir <- "/Users/kristensandys-thomas/Desktop/H3K27ac_diff_peak/test"

input_dir <- file.path(project_dir, "02_top_regions")
output_dir <- file.path(project_dir, "03_great", "raw")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

#defines the upregulated and downregulated files from one comparison.
input_files <- file.path(
  input_dir,
  c(
    "NPM1.2_vs_NPM1.1_NPM1.2_UPREGULATED_regions_top100.tsv",
    "NPM1.2_vs_NPM1.1_NPM1.2_DOWNREGULATED_regions_top100.tsv"
  )
)


# 2. GREAT settings
genome_build <- "hg38"
basal_upstream <- 5000
basal_downstream <- 1000
extension <- 1000000

# 3. Read the genomic peak coordinates
read_regions <- function(input_file) {
  data <- read.delim(
    input_file,
    row.names = 1,
    check.names = FALSE
  )

  # remove white spaces 
  peak_id <- trimws(rownames(data))

  # Split coordinates formatted as "chr:start-end" into three columns.
  chromosome <- sub(":.*", "", peak_id)
  start <- as.integer(sub("-.*", "", sub(".*:", "", peak_id)))
  end <- as.integer(sub(".*-", "", peak_id))

  region_info <- data.frame(
    peak_id,
    chromosome,
    start,
    end
  )

  cbind(region_info, data, row.names = NULL)
}

# 4. Assign genes and calculate distances to the TSS
annotate_regions <- function(regions) {
  # Convert the peak coordinates to a GRanges object.
  peaks <- GRanges(
    seqnames = regions$chromosome,
    ranges = IRanges(regions$start, regions$end)
  )

  # Numeric names allow the rGREAT results to be linked back to input rows.
  names(peaks) <- as.character(seq_len(nrow(regions)))

  # Run GREAT locally using hg38 TSS annotations and the basal-plus-extension
  # regulatory-domain model. GO:BP is required to construct the GreatObject
  great_result <- great(
    gr = peaks,
    gene_sets = "GO:BP",
    tss_source = paste0("GREAT:", genome_build),
    mode = "basalPlusExt",
    extend_from = "TSS",
    basal_upstream = basal_upstream,
    basal_downstream = basal_downstream,
    extension = extension,
    exclude = "gap",
    verbose = FALSE
  )

  # Assign genes and calculate signed distances using each peak midpoint.
  associations <- getRegionGeneAssociations(
    great_result,
    by_middle_points = TRUE,
    use_symbols = TRUE
  )

  # Keep peaks without an associated gene in the output with empty fields.
  regions$associated_genes <- ""
  regions$distance_to_tss_bp <- ""

 # if GREAT found no associations, return the regions with eempty annotations  
  if (length(associations) == 0) return(regions)

  association_rows <- as.integer(names(associations))

  # stop if the GREAT results cannot be linked back to the input rows. 
  if (anyNA(association_rows)) {
    stop("rGREAT did not preserve the numeric peak identifiers.")
  }
  
  # store the genes and the distances 
  genes_by_peak <- as.list(associations$annotated_genes)
  distances_by_peak <- as.list(associations$dist_to_TSS)

  # match each GRAET-associated gene and TSS dist. to the orgional peak
  for (i in seq_along(association_rows)) {
    row <- association_rows[i]
    genes <- as.character(genes_by_peak[[i]])
    distances <- as.integer(distances_by_peak[[i]])

    # Keep each gene once while retaining its paired TSS distance.
    keep <- !duplicated(genes)

    regions$associated_genes[row] <- paste(genes[keep], collapse = ";")
    regions$distance_to_tss_bp[row] <- paste(distances[keep], collapse = ";")
  }

  regions
}

# 5. Process the upregulated and downregulated files
for (input_file in input_files) {
  regions <- read_regions(input_file)
  annotated <- annotate_regions(regions)

  output_name <- sub(
    "\\.tsv$",
    "_GREAT_region_gene_mapping.tsv",
    basename(input_file)
  )

  # store the reuslts. 
  output_file <- file.path(output_dir, output_name)

  write.table(
    annotated,
    output_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )
}
