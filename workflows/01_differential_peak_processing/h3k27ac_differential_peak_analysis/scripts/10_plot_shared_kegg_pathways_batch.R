
# Plot shared KEGG pathways across all five NPM1 subtypes.
# Create separate plots for upregulated and downregulated pathways.

library(ggplot2)

top_n <- 16
top_label <- paste0("top", top_n)

groups <- paste0("NPM1.", 1:5)
directions <- c("UPREGULATED", "DOWNREGULATED")

# 1. Define the production input and output locations
base_dir <- "/Users/kristensandys-thomas/Desktop/H3K27ac_diff_peak/h3k27ac_differential_peak_analysis"

input_dir <- file.path(
  base_dir,
  paste0("compared_", top_label, "_KEGG_pathways_3kb_group_specific")
)

output_dir <- file.path(
  base_dir,
  paste0("KEGG_plots_", top_label, "_3kb_group_specific")
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 2. Plot the upregulated and downregulated shared pathways
for (direction in directions) {
  plot_prefix <- paste0(
    "ALL_NPM_",
    direction,
    "_3kb_GREAT_KEGG_",
    top_label
  )

  input_file <- file.path(
    input_dir,
    paste0(plot_prefix, "_shared_summary.tsv")
  )

  shared <- read.delim(input_file, check.names = FALSE)

  # Convert the five subtype columns into a long table for plotting.
  plot_list <- lapply(groups, function(subtype) {
    data.frame(
      Pathway = shared$Pathway,
      Subtype = subtype,
      Present = shared[[subtype]]
    )
  })

  plot_data <- do.call(rbind, plot_list)

  # Replace YES and NO with clearer plot labels.
  plot_data$Present <- ifelse(
    plot_data$Present == "YES",
    "Present",
    "Absent"
  )

  # Order pathways by the number of subtypes in which they occur.
  pathway_order <- shared$Pathway[
    order(shared$Frequency, shared$Pathway)
  ]

  plot_data$Pathway <- factor(
    plot_data$Pathway,
    levels = pathway_order
  )

  plot_data$Subtype <- factor(
    plot_data$Subtype,
    levels = groups
  )

  # Plot pathway presence and absence across the five NPM1 subtypes.
  plot <- ggplot(
    plot_data,
    aes(x = Subtype, y = Pathway, fill = Present)
  ) +
    geom_tile(color = "white", linewidth = 1) +
    scale_x_discrete(labels = paste0("AML-", groups)) +
    scale_fill_manual(
      values = c(
        "Absent" = "grey92",
        "Present" = "#F4A6C1"
      )
    ) +
    labs(
      title = paste0("Shared ", tolower(direction), " KEGG pathways"),
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

  # Increase the plot height when more shared pathways are present.
  plot_height <- max(8, 2.5 + 0.38 * nrow(shared))

  png_file <- file.path(
    output_dir,
    paste0(plot_prefix, "_presence_absence_plot.png")
  )

  # Write a high-resolution PNG file.
  ggsave(
    png_file,
    plot,
    width = 13,
    height = plot_height,
    dpi = 300,
    limitsize = FALSE
  )
}
