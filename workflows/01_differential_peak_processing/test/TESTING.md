# H3K27ac Differential-Peak Analysis test 

## Purpose

This folder contains the scripts used to test and establish the H3K27ac
differential-peak analysis workflow

## 1. Select the top 100 peaks/ regions from the deseq2 files produced from seq2science ChIP-seq on H3K27ac data. 

peaks == regions 

This stage selects the top 100 peaks in both directions from a DESeq2 result file. 

At this stage, no packages were required, simple base R functions were used. And I ran it from my R workspace on my desktop. 

## Run script: 

``` r
source("01_test_extract_top_100_regions.R")
``` 
I used one input file, called GRCh38-condition_NPM1S2_NPM1S1.diffexp.tsv
I validated the outputs by visually inspecting the top 100 ranked regions selected and compared that to the origional DESeq2 file.
The script correctly selected the ranked 100 regions in both directions. 
A positive log2FC is associated with those regions in NPM1.2 as increased H3K27ac enrichment, and negative log2FC is associated with decreased H3K27ac enrichment in NPM1.2. 
The output file was stored in the folder called `02_top_regions/`. 

## Stage 1 passed. 


## 2. Assign genes to the selected peaks using GREAT

## Input files
Input files: up- and down-regulated top-100 regions from the NPM1.2 vs NPM1.1 comparison generated in Stage 1.
The input files were read from `02_top_regions/`. 

Local rGREAT was used to associate regions with genes. Gene assignments and signed TSS distances were calculated 
from peak midpoints, where negative values indicate upstream regions and positive values indicate downstream regions.

### Run script

```r
setwd("/Users/kristensandys-thomas/Desktop/H3K27ac_diff_peak/test")
source("02_test_annotate_regions_with_great.R")
```

Each input file contained 100 peak observations. Peak identifiers formatted as
`chr:start-end` were separated into chromosome, start, and end coordinates and
converted to `GRanges` objects.

### packages 

rGREAT (version 2.12.2) with `GenomicRanges` and `IRanges`.


### Output 

The annotated peak tables were written to:

```text
03_great/raw/
```

## Validation 
I inspected the files generated. 
To determine if peak co-ordiantes were correctly annotated to genes, I selected a few co-ordinates individually and 
out them into the UCSC genome browser. The peaks were correctly assigned to genes. 

You could also validate the gene annotations by manually submitting the list of genomic co-ordiantes 
to http://great.stanford.edu/public/html/
however this server is temperamental. 

## Stage 2 passed.


## 3. Extract genes and keep the closest TSS association

The upregulated and downregulated GREAT results for NPM1.2 vs NPM1.1 were
processed separately. Peaks without associated genes were removed, and paired
gene–TSS distances were extracted. If a gene was associated with multiple
peaks, only the association with the shortest absolute TSS distance was
retained. Negative distances = indicate upstream and positive distance = downstream
associations. 


### Run script

I first set my working directory to the testing folder, and then ran the script
with source()

```r
setwd("/Users/kristensandys-thomas/Desktop/H3K27ac_diff_peak/test")
source("03_test_extract_genes_with_tss_distance.R")
```

No packages needed, as only base R functions were executed. 

### Input

The Stage 2 GREAT mapping files were read from:

```text
03_great/raw/
```

The two input files contained the upregulated and downregulated top-100 peak
annotations for NPM1.2 vs NPM1.1.

### Output

The cleaned outputs were written to:

```text
03_great/cleaned/
```

Two outputs were created for each direction:
1. A .tsv table containing each unique gene, its signed TSS distance, and its
  absolute TSS distance.
2. A .txt file containing gene names only.

### Validation

The upregulated result contained 147 unique genes, and the downregulated
result contained 128 unique genes. Beucase GREAT can associate one peak to 
multiple gens, more genes than input peaks were observed. 

I checked the distance tables and confirmed that no duplicate genes were
retained. The final tables and gene-only files were ordered alphabetically by
gene name to make manual inspection easier. The gene-only files contained the
same number of genes as their corresponding distance tables.

## Stage 3 passed.


## 4. Select genes within 3 kb of a TSS

This stage tested the selection of GREAT-associated genes located no more than
3,000 bp from a transcription start site. The upregulated and downregulated
NPM1.2 vs NPM1.1 gene-distance tables generated in Stage 3 were processed
separately.

The filter used the absolute TSS distance. A gene was retained when its absolute distance was less or equal 
to 3,000 bp.

### Run script

```r
setwd("~/Desktop/H3K27ac_diff_peak/test")
source("04_test_select_genes_within_3kb.R")
```

### Input

The gene-distance tables were read from:

```text
03_great/cleaned/
```

The upregulated input contained 147 unique genes and the downregulated input
contained 128 unique genes.

### Output

The filtered outputs were written to:

```text
03_great/cleaned_3kb/
```

For each direction, the script produced a gene-distance table and a gene-only
text file. The results were ordered alphabetically by gene name for inspection. 

### Validation
- For NPM1.2, the upregualted genes retained within the 3 Kb were 36 genes, and PRDM16 was incldued as it was 926 
bp away from the TSS.
- NPM1.2's downregulated genes, 98 were retained. 
- All retained genes had an absolute TSS distance no greater than 3,000 bp.
 

## Stage 4 passed.


## 5. Pool candidate genes across comparisons

I tested whether the within-3-kb gene lists could be pooled across
multiple NPM1.2 comparisons. Upregulated and downregulated genes were pooled
separately.

The test included gene lists from GREAT gene associations. 

- NPM1.2 vs NPM1.1
- NPM1.2 vs NPM1.3
- NPM1.2 vs NPM1.4

### Run script

I ran the script from my R worksapce: 

```r
setwd("~/Desktop/H3K27ac_diff_peak/test")
source("06_test_pool_great_genes.R")
```

### Input

The six upregulated and downregulated gene-only files were read from:

```text
03_great/cleaned_3kb/
```

### Output

The pooled gene lists and presence/absence tables were written to:

```text
summary/
```

### Validation

I manually validated the presence/absence and frequency for PRDM16, as I knew it was 
was present in the upregulated NPM1.2 vs NPM1.1 list, with an
absolute TSS distance of 926 bp, and in the upregulated NPM1.2 vs NPM1.4 list.
It was absent from NPM1.2 vs NPM1.3. The pooled reported it's freuqnecy as 2. 

The gene OSTC was present in the downregulated NPM1.2 vs NPM1.1 list, with an
absolute TSS distance of 181 bp, and in the downregulated NPM1.2 vs NPM1.4
list. It was absent from NPM1.2 vs NPM1.3. The pooled reported it's frequencya as 2.. 
Therefore this script was successful. 

## Stage 5 passed.


## 6. Run KEGG enrichment with Enrichr

This stage tested KEGG enrichment of the pooled NPM1.2 candidate genes. The
upregulated and downregulated pooled gene lists from Stage 5 were submitted to
Enrichr separately using the `KEGG_2026` gene-set library.

### Run script

```r
setwd("~/Desktop/H3K27ac_diff_peak/test")
source("06_test_kegg_enrichr.R")
```

The script required the `httr` (version ‘1.4.8’) and `jsonlite`  ( version ‘2.0.0’)
packages and an internet connection.

### Input

The pooled NPM1.2 gene lists were read from:

```text
summary/
```

### Output

The script-generated KEGG enrichment tables were written to:

```text
KEGG/
```

### Manual validation

I manually submitted the same pooled upregulated and downregulated gene lists
to the Enrichr website (https://maayanlab.cloud/Enrichr/enrich) and selected the `KEGG_2026` library. I downloaded the
website results and saved them in testing directory:

```text
KEGG_manual_validation_NPM1.2/
```

I compared the downloaded results with the tables generated by the R script.
I also confirmed that PRDM16 appeared among the overlapping genes in the NPM1.2 upregulated KEGG results. 
It occurred in the THERMOGENESIS and LYSINE DEGRADATION pathways in both the script-generated
and manually downloaded results.

## Stage 6 passed.


## 7. Select the top 16 KEGG pathways

This stage tested the selection of the top 16 pathways from the NPM1.2 KEGG
enrichment results. The upregulated and downregulated enrichment tables were
processed separately. Pathways were ordered by adjusted p-value, followed by
raw p-value when adjusted p-values were equal.

### Run script

```r
setwd("~/Desktop/H3K27ac_diff_peak/test")
source("08_test_filter_top_kegg_pathways.R")
```

No additional packages were required because the script used base R functions.

### Input

The NPM1.2 upregulated and downregulated `KEGG_2026` enrichment tables were
read from:

```text
KEGG/
```

### Output

The selected pathway tables were written to:

```text
KEGG/top16_KEGG_pathways_3kb_group_specific/
```

### Validation

I confirmed that two output files were generated, one for each direction, and
that each table contained the first 16 pathway observations. I also compared the
outputs with the complete Stage 6 KEGG enrichment tables after ordering them by
adjusted p-value and raw p-value. Both outputs contained the expected first 16
pathways in the correct order.

## Stage 7 passed.


## 8. Combine the top 16 KEGG pathways across NPM1 subtypes

This stage tested whether top-16 KEGG pathway tables from different NPM1
subtypes could be combined. The test included NPM1.1, NPM1.2, and NPM1.3, with
upregulated and downregulated results processed separately.

### Run script

```r
setwd("~/Desktop/H3K27ac_diff_peak/test")
source("09_test_combine_top_kegg_pathways.R")
```

No additional packages were required because the script used base R functions.

### Input

The top-16 pathway tables were read from:

```text
KEGG/top16_KEGG_pathways_3kb_group_specific/
```

### Output

The two combined tables were written to:

```text
compared_top16_KEGG_pathways_3kb_group_specific/
```

One table contained the combined upregulated results and the other contained
the combined downregulated results. The script added `NPM_group` and
`Direction` columns so that the source of each pathway row remained clear.

### Validation

I confirmed that each combined table contained 48 rows: 16 pathways for each
of the three test subtypes. The `NPM_group` column contained NPM1.1, NPM1.2,
and NPM1.3 in the expected order. I also confirmed that the `Direction` column
was consistently upregulated or downregulated in the corresponding output.

## Stage 8 passed.


## 9. Filter and identify shared KEGG pathways

This stage tested the filtering of the combined top-16 KEGG tables. Pathways
were retained when they contained at least two overlapping genes. A pathway
was classified as shared when it occurred in at least two of the three test
subtypes: NPM1.1, NPM1.2, and NPM1.3.

### Run script

```r
setwd("~/Desktop/H3K27ac_diff_peak/test")
source("09.1_test_filter_shared_kegg_pathways.R")
```

No additional packages were required because the script used base R functions.

### Input and output

The combined upregulated and downregulated tables were read from and the
filtered results were written to:

```text
compared_top16_KEGG_pathways_3kb_group_specific/
```

The original combined tables were preserved with the suffix
`_combined_with_genes_unfiltered.tsv`. The filtered combined tables included
an `overlapping_gene_count` column. Separate `_shared_summary.tsv` files
recorded pathway presence and frequency across the three test subtypes.

### Validation

I confirmed that every retained pathway contained at least two overlapping
genes. The upregulated table retained 46 pathway rows and the downregulated
table retained 45. The shared summaries each contained seven pathways, and
every shared pathway had a frequency of at least two.

## Stage 9 passed.


## 10. Plot shared KEGG pathway presence across NPM1 subtypes

This stage tested the plotting of shared KEGG pathway presence and absence
across NPM1.1, NPM1.2, and NPM1.3. Upregulated and downregulated pathways were
plotted separately.

### Run script

```r
setwd("~/Desktop/H3K27ac_diff_peak/test")
source("10_test_plot_shared_kegg_pathways.R")
```

The script required the `ggplot2` R package.

### Input

The upregulated and downregulated shared-pathway summaries were read from:

```text
compared_top16_KEGG_pathways_3kb_group_specific/
```

### Output

The pathway presence/absence plots were written to:

```text
KEGG_plots_top16_3kb_group_specific/
```

Pathway presence was displayed in pink and absence in grey. The rows represent
shared pathways and the columns represent the three test NPM1 subtypes.

### Validation

I confirmed that separate upregulated and downregulated PNG files were
generated. 

## Stage 10 passed.


## 11. Plot shared KEGG pathways with genes and significance

This stage tested the addition of overlapping gene labels and adjusted
p-value significance stars to the shared KEGG pathway plots. NPM1.1, NPM1.2,
and NPM1.3 were plotted, with upregulated and downregulated pathways processed
separately.

### Run script

```r
setwd("~/Desktop/H3K27ac_diff_peak/test")
source("11_test_plot_shared_pathways_with_genes.R")
```

The script required the `ggplot2` R package.

### Input

The filtered combined pathway tables were read from:

```text
compared_top16_KEGG_pathways_3kb_group_specific/
```

### Output

The annotated pathway plots were written to:

```text
KEGG_plots_top16_3kb_group_specific/
```

Each present pathway tile displayed up to two overlapping genes. `+N`
indicated the number of additional genes not displayed. Adjusted p-value
significance was represented by `*`, `**`, or `***`.

### Validation

I confirmed that separate upregulated and downregulated PNG files were
generated and that pathway tiles contained the expected gene labels and
significance annotations. Both plots had dimensions of 4,200 by 2,550 pixels.

## Stage 11 passed.


## 12. Trace candidate genes to genomic coordinates

This stage tested whether the pooled NPM1.2 candidate genes could be traced
back to their original GREAT-associated genomic regions. Upregulated and
downregulated genes were processed separately, and only region-gene links
within 3,000 bp of a TSS were retained.

### Run script

```r
setwd("~/Desktop/H3K27ac_diff_peak/test")
source("12_test_trace_genes_to_genomic_coordinates.R")
```

No additional packages were required because the script used base R functions.

### Input

The script used the pooled candidate-gene lists in `summary/`, the raw GREAT
mapping tables in `03_great/raw/`, and the genes retianed from the 3kb threshold in
`03_great/cleaned_3kb/`.

### Output

Separate upregulated and downregulated coordinate tables were written to:

```text
summary/
```

The combined trace table was written to:

```text
TEST_NPM1.2_unique_3kb_GREAT_candidate_genes_with_genomic_coordinates.tsv
```

### Validation
Every retained association was within 3,000 bp of a
TSS. I manually checked a few genes like PRDM16 in the results with the Raw GREAT results,
and the correct genomic co-ordiantes as well, p-values and gene distances were correct. 
The same went for various other genes I tested. 

## Stage 12 passed.
