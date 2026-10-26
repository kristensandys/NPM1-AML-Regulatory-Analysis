# Plot shared KEGG pathways with overlapping genes and significance.
# Use the NPM1.1, NPM1.2, and NPM1.3 test results.

library(ggplot2)

top_n <- 16
top_label <- paste0("top", top_n)
min_groups <- 2

subtypes <- c("NPM1.1", "NPM1.2", "NPM1.3")
directions <- c("UPREGULATED", "DOWNREGULATED")

present_colour <- "#F4A6C1"
absent_colour <- "grey92"

# 1. Define the input and output locations
project_dir <- "/Users/kristensandys-thomas/Desktop/H3K27ac_diff_peak/test"

input_dir <- file.path(
  project_dir,
  paste0("compared_", top_label, "_KEGG_pathways_3kb_group_specific")
)

output_dir <- file.path(
  project_dir,
  paste0("KEGG_plots_", top_label, "_3kb_group_specific")
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 2. Separate the semicolon-delimited overlapping genes
split_genes <- function(x) {
  if (is.na(x) || trimws(x) == "") return(character())

  genes <- trimws(unlist(strsplit(as.character(x), ";", fixed = TRUE)))
  unique(genes[genes != "" & genes != "NA"])
}

# 3. Show up to two genes in each pathway tile
make_gene_label <- function(x) {
  genes <- split_genes(x)
  if (length(genes) == 0) return("")

  shown <- head(genes, 2)
  label <- paste(shown, collapse = "\n")
  extra <- length(genes) - length(shown)

  if (extra > 0) label <- paste0(label, "\n+", extra)

  label
}

# 4. Add significance stars from the adjusted p-value
make_stars <- function(padj) {
  if (is.na(padj)) return("")
  if (padj < 0.001) return("***")
  if (padj < 0.01) return("**")
  if (padj < 0.05) return("*")
  ""
}

# 5. Plot the upregulated and downregulated pathways separately
for (direction in directions) {
  plot_prefix <- paste0(
    "TEST_ALL_NPM_",
    direction,
    "_3kb_GREAT_KEGG_",
    top_label
  )

  input_file <- file.path(
    input_dir,
    paste0(plot_prefix, "_combined_with_genes.tsv")
  )

  combined <- read.delim(input_file, check.names = FALSE)

  pathway_data <- data.frame(
    Pathway = combined$term,
    Subtype = combined$NPM_group,
    Genes = combined$overlapping_genes,
    Gene_count = combined$overlapping_gene_count,
    Adjusted_p_value = as.numeric(combined$adjusted_p_value)
  )

  # Keep pathways found in at least two test subtypes.
  pathway_subtypes <- unique(pathway_data[c("Pathway", "Subtype")])
  pathway_frequency <- table(pathway_subtypes$Pathway)

  shared_pathways <- names(
    pathway_frequency[pathway_frequency >= min_groups]
  )

  pathway_data <- pathway_data[
    pathway_data$Pathway %in% shared_pathways,
  ]

  # Add rows for subtype-pathway combinations that are absent.
  plot_data <- expand.grid(
    Pathway = shared_pathways,
    Subtype = subtypes
  )

  plot_data <- merge(
    plot_data,
    pathway_data,
    by = c("Pathway", "Subtype"),
    all.x = TRUE,
    sort = FALSE
  )

  plot_data$Gene_count[is.na(plot_data$Gene_count)] <- 0
  plot_data$Genes[is.na(plot_data$Genes)] <- ""

  plot_data$Present <- ifelse(
    plot_data$Gene_count > 0,
    "Present",
    "Absent"
  )

  plot_data$Present <- factor(
    plot_data$Present,
    levels = c("Absent", "Present")
  )

  plot_data$Frequency <- as.integer(
    pathway_frequency[plot_data$Pathway]
  )

  # Put PRDM16 first in the selected thermogenesis tiles when present.
  highlight_rows <- direction == highlight_direction &
    toupper(plot_data$Pathway) == highlight_pathway &
    plot_data$Subtype %in% highlight_subtypes

  for (i in which(highlight_rows)) {
    genes <- split_genes(plot_data$Genes[i])

    if (highlight_gene %in% genes) {
      genes <- c(highlight_gene, genes[genes != highlight_gene])
      plot_data$Genes[i] <- paste(genes, collapse = ";")
    }
  }

  plot_data$Label <- vapply(
    plot_data$Genes,
    make_gene_label,
    character(1)
  )

  significance <- vapply(
    plot_data$Adjusted_p_value,
    make_stars,
    character(1)
  )

  has_stars <- significance != ""

  plot_data$Label[has_stars] <- paste0(
    plot_data$Label[has_stars],
    "\n",
    significance[has_stars]
  )

  plot_data$Label[plot_data$Present == "Absent"] <- ""

  # Order pathways by how many subtypes contain them.
  pathway_order <- unique(plot_data[c("Pathway", "Frequency")])
  pathway_order <- pathway_order[
    order(pathway_order$Frequency, pathway_order$Pathway),
  ]

  plot_data$Pathway <- factor(
    plot_data$Pathway,
    levels = pathway_order$Pathway
  )

  plot_data$Subtype <- factor(
    plot_data$Subtype,
    levels = subtypes
  )

  plot <- ggplot(
    plot_data,
    aes(x = Subtype, y = Pathway, fill = Present)
  ) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(
      aes(label = Label),
      size = 3,
      fontface = "bold",
      lineheight = 0.9,
      color = "black"
    ) +
    scale_x_discrete(labels = paste0("AML-", subtypes)) +
    scale_fill_manual(
      values = c(
        "Absent" = absent_colour,
        "Present" = present_colour
      ),
      drop = FALSE
    ) +
    labs(
      title = paste0("Shared ", tolower(direction), " KEGG pathways"),
      caption = paste0(
        "+N indicates additional genes; ",
        "* padj < 0.05, ** padj < 0.01, *** padj < 0.001"
      ),
      x = NULL,
      y = NULL,
      fill = NULL
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 20),
      axis.text.x = element_text(face = "bold", size = 13),
      axis.text.y = element_text(size = 11),
      panel.grid = element_blank(),
      legend.position = "bottom"
    )

  plot_height <- max(8.5, 2.5 + 0.65 * length(shared_pathways))

  output_file <- file.path(
    output_dir,
    paste0(plot_prefix, "_shared_pathways_with_genes.png")
  )

  ggsave(
    output_file,
    plot,
    width = 14,
    height = plot_height,
    dpi = 300,
    limitsize = FALSE
  )
}
