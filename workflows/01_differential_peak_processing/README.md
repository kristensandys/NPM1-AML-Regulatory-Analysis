# H3K27ac differential-peak analysis

Research scripts used during my MSc internship to analyse differential H3K27ac ChIP-seq peaks across five NPM1 groups.

## Contents

- [Workflow documentation](h3k27ac_differential_peak_analysis/README.md): analysis stages and the original execution record.
- [Production scripts](h3k27ac_differential_peak_analysis/scripts): peak selection, local rGREAT annotation, candidate-gene filtering within 3 kb of a TSS, pooling, KEGG enrichment, pathway comparisons, plotting, and coordinate tracing.
- [Validation documentation](test/TESTING.md) and validation scripts in `test/`.
- [Example input](examples/): one differential-peak input table.
- [Example figure](figures/): a representative analysis plot.

## Scope of this code release

This repository includes code, documentation, one example input table, and one representative figure. Other research inputs, metadata, generated results, archived script versions,
and R session history are excluded. References to existing results and passed checks in the original documentation describe the author's prior runs; they do not mean those results are bundled here or that this copy has been rerun.

## Before running

The scripts retain the paths used in the original analysis. Update `base_dir` in the R scripts and `new_root` / `old_root` in the shell scripts for your environment. Validation scripts also require their original paths and input fixtures to be adapted. Restore the required differential-peak inputs to the expected comparison directories before running the workflow.

The workflow starts from existing DESeq2 differential-peak results; it does not perform raw-read processing or differential testing. Peaks are selected by log2 fold change without an adjusted-p-value selection threshold. The `UPREGULATED` and `DOWNREGULATED` labels describe ChIP-seq signal direction and do not alone establish gene-expression changes.

Required R packages include `rGREAT`, `GenomicRanges`, `IRanges`, `httr`, `jsonlite`, and `ggplot2`. The Enrichr stage requires internet access and the configured KEGG library. Shell setup scripts use Bash associative arrays and therefore require Bash 4 or later.

The NPM1.5-versus-NPM1.4 differential result is absent from the original workflow, as recorded in the workflow documentation.

