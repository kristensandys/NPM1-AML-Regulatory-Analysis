# Run KEGG enrichment on the pooled NPM1.2 candidate gene lists whe NPM1.2 was compared 
# with NPM1.1, NPM1.3, and NPM1.4. 


library(httr)
library(jsonlite)

# 1. project location 
project_dir <- "/Users/kristensandys-thomas/Desktop/H3K27ac_diff_peak/test"

input_dir <- file.path(project_dir, "summary")
output_dir <- file.path(project_dir, "KEGG")

# creates the new directory 
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

input_files <- file.path(
  input_dir,
  c(
    "NPM1.2_UPREGULATED_unique_3kb_GREAT_candidate_genes.txt",
    "NPM1.2_DOWNREGULATED_unique_3kb_GREAT_candidate_genes.txt"
  )
)

# 2. Define the Enrichr settings
enrichr_url <- "https://maayanlab.cloud/Enrichr"
kegg_library <- "KEGG_2026" # gene library used 

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
 
   # submit the genes as newline-separated list. 
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
  # stop if the submission failed. 
  stop_for_status(response)
  # wait before sending next request.
  Sys.sleep(pause_seconds)
  
  # convert Enrichr's response into an R object.  
  fromJSON(content(response, as = "text", encoding = "UTF-8"))
}

# 5. Retrieve the KEGG enrichment results from Enrichr
get_kegg <- function(user_list_id) {
  
  # Request results for the submitted gene list & KEGG library
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
  # stop if the request failed, then pause before the next request. 
  stop_for_status(response)
  Sys.sleep(pause_seconds)

  # convert the returned JSON into an R object 
  results <- fromJSON(
    content(response, as = "text", encoding = "UTF-8"),
    simplifyVector = FALSE
  )

  # check existed result exists
  if (is.null(results[[kegg_library]])) {
    stop("No results were returned for ", kegg_library)
  }
  # return the results for the selected KEGG library. 
  results[[kegg_library]]
}

# 6. Convert the Enrichr results into a table
format_results <- function(enrichment) {
  
  # return an empty table with the correct columns if no pathways were found. 
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

  # convert enriched pathways into one row for each 
  rows <- lapply(enrichment, function(x) {
    data.frame(
      rank = as.integer(x[[1]]),
      term = as.character(x[[2]]),
      p_value = as.numeric(x[[3]]),
      odds_ratio = as.numeric(x[[4]]),
      combined_score = as.numeric(x[[5]]),
      
      # store overlapping genes into one field, separated by ";" 
      overlapping_genes = paste(unlist(x[[6]]), collapse = ";"),
      adjusted_p_value = as.numeric(x[[7]])
    )
  })
  
  #combine resylts into one table 
  do.call(rbind, rows)
}

# 7. Run enrichment for the upregulated and downregulated gene lists
for (input_file in input_files) {
  # get the file without the .txt extension
  description <- tools::file_path_sans_ext(basename(input_file))
  prefix <- paste0(description, "_", kegg_library)

  # define name and location of the output file. 
  results_file <- file.path(
    output_dir,
    paste0(prefix, "_enrichment.tsv")
  )

  # read the genes from the current pooled gene list. 
  genes <- read_genes(input_file)

  # submit genes into Enrichr & save the submission ID. 
  submission <- submit_genes(genes, description)
  user_list_id <- submission$userListId

  enrichment <- get_kegg(user_list_id)
  results <- format_results(enrichment)

  #write the KEGG enrichment reuslts to a .tsv
  write.table(
    results,
    results_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )
}
