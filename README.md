
<!-- README.md is generated from README.Rmd. Please edit that file -->

# markeR

<!-- badges: start -->

![Development
Status](https://img.shields.io/badge/status-development-yellowgreen)
<!-- badges: end -->

**markeR** provides a suite of methods for using gene sets (signatures)
to quantify and evaluate the extent to which a given gene signature
marks a specific phenotype. The package implements various scoring,
enrichment and classification approaches, along with tools to compute
performance metrics and visualize results, making it a valuable resource
for transcriptomics research (bulk RNA-seq).

## Installation

You can install the development version of markeR from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("DiseaseTranscriptomicsLab/markeR")
```

## Methods

The package implements three main approaches for gene signature
analysis:

**1. Score-based Methods:**

These approaches compute a numerical score representing the expression
of a gene signature per sample.

-   **Log2 Median-Centered Scoring:** Measures signature expression by
    centering expression values relative to their median.
-   **Single-Sample Gene Set Enrichment Analysis (ssGSEA):** Calculates
    an enrichment score for a given gene set in each sample.

**2. Enrichment-based Methods:**

These methods evaluate whether a gene signature is significantly
enriched in a ranked gene list.

-   **Gene Set Enrichment Analysis (GSEA):** A rank-based method for
    testing enrichment significance in predefined gene sets.

**3. Classification-based Methods:**

These approaches classify samples into phenotypic groups based on gene
signature expression.

-   **Random Forest Classification:** Uses decision trees to predict
    sample labels based on gene signature scores.

## Main Functions and Future Modules

The current release of **markeR** includes two primary functions for
score-based analysis:

-   **CalculateScores:** Calculates gene signature scores for each
    sample using either the ssGSEA or log2 median-centered method.
-   **PlotScores:** Visualizes the calculated scores across conditions
    using violin plots.

Future updates will expand the package with additional pairs of
functions to:

-   Perform enrichment-based analysis (to calculate and visualize GSEA
    statistics).
-   Conduct classification-based analysis (to train classifiers, e.g.,
    using Random Forests, and to evaluate performance via ROC curves).
-   Provide gene-level visualization modules (e.g., heatmaps) to display
    expression patterns of individual genes within signatures.

## Example

This example demonstrates the calculation of a log2-median-centered
score using mock RNA-seq expression data. The dataset is derived from
the Marthandan et al. (2016) study (GSE63577) and includes fibroblast
samples under replicative senescent and proliferative conditions. (see
`?counts_example` and `?metadata_example` for more details).This example
showcases how to compute gene expression scores by applying the
log2-median-centered transformation.

We are using as an example a very simple senescence gene signature,
composed of the genes usually associated with senescence.

``` r
library(markeR)
```

``` r
# Define simple Senescence Signature
SimpleSenescenceSignature <- c("CDKN1A", "CDKN2A", "GLB1","TP53","CCL2")
```

``` r
data(metadata_example)
data(counts_example)

# Load example data
head(metadata_example)
#>       sampleID      DatasetID   CellType     Condition       SenescentType
#> 252 SRR1660534 Marthandan2016 Fibroblast     Senescent Telomere shortening
#> 253 SRR1660535 Marthandan2016 Fibroblast     Senescent Telomere shortening
#> 254 SRR1660536 Marthandan2016 Fibroblast     Senescent Telomere shortening
#> 255 SRR1660537 Marthandan2016 Fibroblast Proliferative                none
#> 256 SRR1660538 Marthandan2016 Fibroblast Proliferative                none
#> 257 SRR1660539 Marthandan2016 Fibroblast Proliferative                none
#>                         Treatment
#> 252 PD72 (Replicative senescence)
#> 253 PD72 (Replicative senescence)
#> 254 PD72 (Replicative senescence)
#> 255                         young
#> 256                         young
#> 257                         young
counts_example[1:5,1:5]
#>          SRR1660534 SRR1660535 SRR1660536 SRR1660537 SRR1660538
#> A1BG        9.94566   9.476768   8.229231   8.515083   7.806479
#> A1BG-AS1   12.08655  11.550303  12.283976   7.580694   7.312666
#> A2M        77.50289  56.612839  58.860268   8.997624   6.981857
#> A4GALT     14.74183  15.226083  14.815891  14.675780  15.222488
#> AAAS       47.92755  46.292377  43.965972  47.109493  47.213739
```

``` r
df_Scores <- CalculateScores(data = counts_example, 
                             metadata = metadata_example, 
                             method = "logmedian", 
                             gene_sets = list(Senescence=SimpleSenescenceSignature))

senescence_triggers_colors <- c(
  "none" = "#E57373",  # Soft red   
  "Telomere shortening" = "#4FC3F7"  # Vivid sky blue  
)
 
plt <- PlotScores(ResultsList = df_Scores, 
                  ColorVariable = "SenescentType", 
                  GroupingVariable="Condition",  
                  method ="logmedian", 
                  ColorValues = senescence_triggers_colors, 
                  ConnectGroups=TRUE, 
                  ncol = NULL, 
                  nrow = NULL, 
                  widthTitle=20, 
                  y_limits = NULL, 
                  legend_nrow = 1, 
                  pointSize=4)  

ggpubr::annotate_figure(plt, top = grid::textGrob("Marthandan et al. 2016", gp = grid::gpar(cex = 1.3, fontsize = 12)))
```

<img src="man/figures/README-exampleScore-1.png" width="40%" />
