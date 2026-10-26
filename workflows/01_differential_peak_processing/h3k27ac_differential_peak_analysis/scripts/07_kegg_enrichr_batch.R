# Run KEGG enrichment on the pooled candidate gene lists from all five NPM1
# subtype projects. Process upregulated and downregulated genes separately.

library(httr)
library(jsonlite)

# 1. the five project locations
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

# 2. Define the Enrichr settings
enrichr_url <- "https://maayanlab.cloud/Enrichr"
kegg_library <- "KEGG_2026" #gene library used 

# Wait between requests and retry if the server is temporarily unavailable.
pause_seconds <- 5
max_attempts <- 8

# 3. Read one pooled gene list and remove whitespace, blank lines, and duplicates.
read_genes <- function(input_file) {
  genes <- trimws(readLines(input_file))
  unique(genes[genes != ""])
}

# 4. Submit one gene list to Enrichr and return its submission ID
submit_genes <- function(genes, description) {
  response <- RETRY(
    "POST",
    paste0(enrichr_url, "/addList"),
    body = list(
      list = paste(genes, collapse = "\n"),
      description = description
    ),
    encode = "multipart",
    times = max_attempts,
    pause_base = pause_seconds,
    pause_cap = 120,
    quiet = TRUE
  )

  stop_for_status(response)
  Sys.sleep(pause_seconds)

  fromJSON(content(response, as = "text", encoding = "UTF-8"))
}

# 5. Retrieve the KEGG enrichment results from Enrichr
get_kegg <- function(user_list_id) {
  response <- RETRY(
    "GET",
    paste0(enrichr_url, "/enrich"),
    query = list(
      userListId = user_list_id,
      backgroundType = kegg_library
    ),
    times = max_attempts,
    pause_base = pause_seconds,
    pause_cap = 120,
    quiet = TRUE
  )

  stop_for_status(response)
  Sys.sleep(pause_seconds)

  results <- fromJSON(
    content(response, as = "text", encoding = "UTF-8"),
    simplifyVector = FALSE
  )

  if (is.null(results[[kegg_library]])) {
    stop("No results were returned for ", kegg_library)
  }

  results[[kegg_library]]
}

# 6. Convert the Enrichr results into a table
format_results <- function(enrichment) {
  if (length(enrichment) == 0) {
    return(
      data.frame(
        rank = integer(),
        term = character(),
        p_value = numeric(),
        odds_ratio = numeric(),
        combined_score = numeric(),
        overlapping_genes = character(),
        adjusted_p_value = numeric()
      )
    )
  }

  rows <- lapply(enrichment, function(x) {
    data.frame(
      rank = as.integer(x[[1]]),
      term = as.character(x[[2]]),
      p_value = as.numeric(x[[3]]),
      odds_ratio = as.numeric(x[[4]]),
      combined_score = as.numeric(x[[5]]),
      overlapping_genes = paste(unlist(x[[6]]), collapse = ";"),
      adjusted_p_value = as.numeric(x[[7]])
    )
  })

  do.call(rbind, rows)
}

# 7. Process the two pooled gene lists from every NPM1 project
for (subtype in names(project_dirs)) {
  project_dir <- project_dirs[[subtype]]

  input_dir <- file.path(project_dir, "summary")
  output_dir <- file.path(project_dir, "KEGG")

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  input_files <- file.path(
    input_dir,
    c(
      paste0(
        subtype,
        "_UPREGULATED_unique_3kb_GREAT_candidate_genes.txt"
      ),
      paste0(
        subtype,
        "_DOWNREGULATED_unique_3kb_GREAT_candidate_genes.txt"
      )
    )
  )

  for (input_file in input_files) {
    description <- tools::file_path_sans_ext(basename(input_file))
    prefix <- paste0(description, "_", kegg_library)

    results_file <- file.path(
      output_dir,
      paste0(prefix, "_enrichment.tsv")
    )

    genes <- read_genes(input_file)

    submission <- submit_genes(genes, description)
    user_list_id <- submission$userListId

    enrichment <- get_kegg(user_list_id)
    results <- format_results(enrichment)

    write.table(
      results,
      results_file,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      na = "NA"
    )
  }
}
