
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

## Table of Contents

-   [Installation](#installation)
-   [Main Functions and Future
    Modules](#main-functions-and-future-modules)
-   [Example](#example)
    -   [Visualise Individual Genes from Senescence
        Signature](#visualise-individual-genes-from-senescence-signature)
        -   [Expression Heatmap](#expression-heatmap)
        -   [Expression Violins](#expression-violins)
        -   [Correlation Heatmap](#correlation-heatmap)
        -   [ROC and AUC](#roc-and-auc)
        -   [Cohen’s D](#cohens-d)
        -   [PCA with only genes of
            interest](#pca-with-only-genes-of-interest)
    -   [Calculate Senescence Scores](#calculate-senescence-scores)
        -   [logmedian method](#logmedian-method)
        -   [ssGSEA method](#ssgsea-method)
        -   [Ranking method](#ranking-method)
        -   [All methods](#all-methods)

## Installation

You can install the development version of markeR from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("DiseaseTranscriptomicsLab/markeR")
```

<!-- ## Methods -->
<!-- The package implements three main approaches for gene signature analysis: -->
<!-- **1. Score-based Methods:**   -->
<!-- These approaches compute a numerical score representing the expression of a gene signature per sample. -->
<!-- - **Log2 Median-Centered Scoring:** Measures signature expression by centering expression values relative to their median. -->
<!-- - **Single-Sample Gene Set Enrichment Analysis (ssGSEA):** Calculates an enrichment score for a given gene set in each sample. -->
<!-- **2. Enrichment-based Methods:**   -->
<!-- These methods evaluate whether a gene signature is significantly enriched in a ranked gene list. -->
<!-- - **Gene Set Enrichment Analysis (GSEA):** A rank-based method for testing enrichment significance in predefined gene sets. -->
<!-- **3. Classification-based Methods:**   -->
<!-- These approaches classify samples into phenotypic groups based on gene signature expression. -->
<!-- - **Random Forest Classification:** Uses decision trees to predict sample labels based on gene signature scores. -->

## Main Functions and Future Modules

The current release of **markeR** includes two primary functions for
score-based analysis:

-   **CalculateScores:** Calculates gene signature scores for each
    sample using either the ssGSEA, log2 median-centered or ranking
    method.
-   **PlotScores:** Calculates and displays the calculated scores across
    conditions using violin plots, density plots or heatmaps, depending
    on the chosen parameters.

It also includes some functions for visualising individual genes from a
gene signature:

-   **IndividualGenes\_Violins:** creates violin plots of gene
    expression data with jittered points and optional faceting, allowing
    for visualization of individual gene expression distributions across
    sample groups.
-   **CorrelationHeatmap:** computes and visualizes a correlation
    heatmap for a given set of genes. Optionally, the heatmap can be
    generated separately for different conditions based on metadata.
-   **ExpressionHeatmap:** generates an expression heatmap with
    customizable sample annotations for a given set of genes.
-   **ROCandAUCplot:** computes ROC curves and AUC values for each gene
    based on gene expression data and sample metadata. It can generate
    ROC plots, an AUC heatmap, or both arranged side‐by‐side.
-   **CohenDHeatmap:** computes Cohen’s d for each gene based on gene
    expression data and sample metadata. The resulting effect sizes are
    then visualized as a heatmap.
-   **plotPCA:** performs PCA on a given dataset and visualizes the
    results using ggplot2. It allows users to specify genes of interest
    (to understand if they are sufficient to explain the main variance
    in the data), customize scaling and centering, and color points
    based on a metadata variable.

Future updates will expand the package with additional pairs of
functions to:

-   Perform enrichment-based analysis (to calculate and visualize GSEA
    statistics).
-   Conduct classification-based analysis (to train classifiers, e.g.,
    using Random Forests, and to evaluate performance via ROC curves).
-   Provide additional gene-level visualization modules to display
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
SimpleSenescenceSignature <- c("CDKN1A", "CDKN2A", "GLB1","TP53","CCL2", "LMNB1", "MKI67" )
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

### Visualise Individual Genes from Senescence Signature

#### Expression Heatmap

``` r
annotation_colors <- list( 
  Condition = c(
    "Senescent"     = "#65AC7C",  # Example color: greenish
    "Proliferative" = "#5F90D4"  # Example color: blueish 
  )
)

ExpressionHeatmap(data=counts_example, 
                  metadata = metadata_example, 
                  genes=SimpleSenescenceSignature,  
                  annotate.by = c("Condition"),
                  annotation_colors = annotation_colors,
                  colorlist = list(low = "#3F4193", mid = "#F9F4AE", high = "#B44141"),
                  cluster_rows = TRUE, 
                  cluster_columns = FALSE,
                  title = "Senescence Genes", 
                  titlesize = 20,
                  legend_position = "right",
                  scale_position="right")
```

<img src="man/figures/README-example_exprheatmap-1.png" width="70%" />

#### Expression Violins

``` r
senescence_triggers_colors <- c(
  "none" = "#E57373",  # Soft red   
  "Telomere shortening" = "#4FC3F7"  # Vivid sky blue  
)


IndividualGenes_Violins(data = counts_example, 
                        metadata = metadata_example, 
                        genes = SimpleSenescenceSignature, 
                        GroupingVariable = "Condition", 
                        plot=T, 
                        ncol=NULL, 
                        nrow=1, 
                        divide=NULL, 
                        invert_divide=FALSE,
                        ColorValues=senescence_triggers_colors, 
                        pointSize=2, 
                        ColorVariable="SenescentType", 
                        title="Senescence Genes", 
                        widthTitle=16,
                        y_limits = NULL,
                        legend_nrow=1, 
                        xlab="Condition",
                        colorlab="") 
```

<img src="man/figures/README-exampleviolins-1.png" width="100%" />

#### Correlation Heatmap

``` r
CorrelationHeatmap(data=counts_example, 
                   metadata = metadata_example, 
                   genes=SimpleSenescenceSignature, 
                   separate.by = "Condition", 
                   method = "spearman",  
                   colorlist = list(low = "#3F4193", mid = "#F9F4AE", high = "#B44141"),
                   limits_colorscale = c(-1,0,1), 
                   widthTitle = 16, 
                   title = "Senescence Genes", 
                   cluster_rows = TRUE, 
                   cluster_columns = TRUE,  
                   detailedresults = FALSE, 
                   legend_position="right",
                   titlesize=20)
```

<img src="man/figures/README-example_heatmap-1.png" width="70%" />

#### ROC and AUC

``` r
senescence_triggers_colors <- c(
  "none" = "#E57373",  # Soft red   
  "Telomere shortening" = "#4FC3F7"  # Vivid sky blue  
)
 
ROCandAUCplot(counts_example, 
              metadata_example, 
              condition_var = "Condition", 
              class = "Senescent", 
              group_var=NULL,
              genes=SimpleSenescenceSignature, 
              plot_type = "all",
              heatmap_params = list(col = list( "#F9F4AE" ,"#B44141"),
                                    limits = c(0.5,1),
                                    cluster_rows=T),
              roc_params = list(nrow=3,
                                ncol=3,
                                colors=senescence_triggers_colors),
              commomplot_params = list(widths=c(0.5,0.3)))
```

<img src="man/figures/README-rocAUCexample-1.png" width="100%" />

#### Cohen’s D

``` r
CohenD_IndividualGenes(counts_example, 
              metadata_example, 
              genes=SimpleSenescenceSignature,
              condition_var = "Condition", 
              class = "Senescent", 
              group_var = NULL,  
              heatmap_params = list(col = list( "#F9F4AE" ,"#B44141"),
                                    limits = NULL,
                                    cluster_rows=T))
```

<img src="man/figures/README-cohendexample-1.png" width="70%" />

#### PCA with only genes of interest

``` r
annotation_colors <- c(  
    "Senescent"     = "#65AC7C",  # Example color: greenish
    "Proliferative" = "#5F90D4"  # Example color: blueish 
)

 
plotPCA(data = counts_example, 
        metadata = metadata_example, 
        genes=SimpleSenescenceSignature, 
        scale=FALSE, 
        center=TRUE, 
        PCs=list(c(1,2), c(2,3), c(3,4)), 
        ColorVariable="Condition",
        ColorValues=annotation_colors,
        pointSize=5,
        legend_nrow=1, 
        ncol=3, 
        nrow=NULL)
```

<img src="man/figures/README-pca-1.png" width="90%" />

### Calculate Senescence Scores

#### logmedian method

The following example uses the **“logmedian”** method for score
calculation.

The user can chose to calculate the gene signature score for each sample
based on one or more predefined gene sets (signatures). If a single
method is chosen, a data frame containing the calculated scores for each
gene signature, including metadata if provided. If method = “all” (see
below for an example), a list is returned where each element corresponds
to a scoring method and contains the respective data frame of scores.

``` r
df_Scores <- CalculateScores(data = counts_example,
                             metadata = metadata_example,
                             method = "logmedian",
                             gene_sets = list(Senescence=SimpleSenescenceSignature))
#> Considering unidirectional gene signature mode for signature Senescence

head(df_Scores$Senescence)
#>       sample      score      DatasetID   CellType     Condition
#> 1 SRR1660534 -0.6894748 Marthandan2016 Fibroblast     Senescent
#> 2 SRR1660535 -0.4483299 Marthandan2016 Fibroblast     Senescent
#> 3 SRR1660536 -0.4596502 Marthandan2016 Fibroblast     Senescent
#> 4 SRR1660537 -0.2198753 Marthandan2016 Fibroblast Proliferative
#> 5 SRR1660538 -0.2672930 Marthandan2016 Fibroblast Proliferative
#> 6 SRR1660539 -0.2623188 Marthandan2016 Fibroblast Proliferative
#>         SenescentType                     Treatment
#> 1 Telomere shortening PD72 (Replicative senescence)
#> 2 Telomere shortening PD72 (Replicative senescence)
#> 3 Telomere shortening PD72 (Replicative senescence)
#> 4                none                         young
#> 5                none                         young
#> 6                none                         young
```

The user can also chose to directly plot the scores.

``` r
# df_Scores <- CalculateScores(data = counts_example, 
#                              metadata = metadata_example, 
#                              method = "logmedian", 
#                              gene_sets = list(Senescence=SimpleSenescenceSignature))

senescence_triggers_colors <- c(
  "none" = "#E57373",  # Soft red   
  "Telomere shortening" = "#4FC3F7"  # Vivid sky blue  
)

cond_cohend <- list(A=c("Senescent"),  
                    B=c("Proliferative"))

PlotScores(data = counts_example, 
           metadata = metadata_example, 
           gene_sets = list(Senescence=SimpleSenescenceSignature),
           ColorVariable = "SenescentType", 
           GroupingVariable="Condition",  
           method ="logmedian", 
           ColorValues = senescence_triggers_colors, 
           ConnectGroups=TRUE, 
           ncol = NULL, 
           nrow = NULL, 
           widthTitle=24, 
           limits = NULL, 
           legend_nrow = 1, 
           pointSize=4,
           cond_cohend=cond_cohend,
           title="Marthandan et al. 2016",
           labsize=7, 
           titlesize = 10)  
```

<img src="man/figures/README-exampleScore-1.png" width="40%" />

Given that some of the genes are expected to be upregulated, while
others to be downregulated in senescence, we can also consider a
*bidirectional signature*.

``` r
SimpleSenescenceSignature_bidirectional <- data.frame(gene=c("CDKN1A", "CDKN2A", "GLB1","TP53","CCL2", "LMNB1", "MKI67" ),
                                                      enrichment=c(1,1,1,1,1,-1,-1))
# 
# df_Scores <- CalculateScores(data = counts_example, 
#                              metadata = metadata_example, 
#                              method = "logmedian", 
#                              gene_sets = list(Senescence=SimpleSenescenceSignature_bidirectional))

PlotScores(data = counts_example, 
           metadata = metadata_example, 
           gene_sets = list(Senescence=SimpleSenescenceSignature_bidirectional),
           ColorVariable = "SenescentType", 
           GroupingVariable="Condition",  
           method ="logmedian", 
           ColorValues = senescence_triggers_colors, 
           ConnectGroups=TRUE, 
           ncol = NULL, 
           nrow = NULL, 
           widthTitle=24, 
           limits = NULL, 
           legend_nrow = 1, 
           pointSize=4,
           cond_cohend=cond_cohend,
           title="Marthandan et al. 2016",
           labsize=7, 
           titlesize = 10)  
```

<img src="man/figures/README-exampleScore_bidirectional-1.png" width="40%" />

For users interested in viewing the overall distribution of scores,
simply omit the `GroupingVariable` or `metadata` parameters. In this
case, the function will automatically generate a grid of density plots,
with each gene signature represented by its own plot.

``` r
PlotScores(data = counts_example, 
           metadata = metadata_example, 
           gene_sets = list(Senescence_Bidirectional = SimpleSenescenceSignature_bidirectional,
                          Senescence  = SimpleSenescenceSignature), 
           method ="logmedian", 
           ColorValues = NULL,  
           ncol = NULL, 
           nrow = NULL, 
           widthTitle=24, 
           limits = NULL,  
           title="Marthandan et al. 2016",
           labsize=8, 
           titlesize = 10)  
```

<img src="man/figures/README-plotscores_density-1.png" width="80%" />

#### ssGSEA method

The following example uses the **“ssGSEA”** method for score
calculation, both for unidirectional and bidirectional signatures.

``` r
#  
# df_Scores <- CalculateScores(data = counts_example, 
#                              metadata = metadata_example, 
#                              method = "ssGSEA", 
#                              gene_sets = list(Senescence=SimpleSenescenceSignature))

senescence_triggers_colors <- c(
  "none" = "#E57373",  # Soft red   
  "Telomere shortening" = "#4FC3F7"  # Vivid sky blue  
)

cond_cohend <- list(A=c("Senescent"),  
                    B=c("Proliferative"))

PlotScores(data = counts_example, 
           metadata = metadata_example, 
           gene_sets = list(Senescence=SimpleSenescenceSignature),
           ColorVariable = "SenescentType", 
           GroupingVariable="Condition",  
           method ="ssGSEA", 
           ColorValues = senescence_triggers_colors, 
           ConnectGroups=TRUE, 
           ncol = NULL, 
           nrow = NULL, 
           widthTitle=24, 
           limits = NULL, 
           legend_nrow = 1, 
           pointSize=4,
           cond_cohend=cond_cohend,
           title="Marthandan et al. 2016",
           labsize=7, 
           titlesize = 10)  
```

<img src="man/figures/README-exampleScoresGSEA_uni-1.png" width="40%" />

``` r
 
# df_Scores <- CalculateScores(data = counts_example, 
#                              metadata = metadata_example, 
#                              method = "ssGSEA", 
#                              gene_sets = list(Senescence=SimpleSenescenceSignature_bidirectional))

senescence_triggers_colors <- c(
  "none" = "#E57373",  # Soft red   
  "Telomere shortening" = "#4FC3F7"  # Vivid sky blue  
)

cond_cohend <- list(A=c("Senescent"),  
                    B=c("Proliferative"))

PlotScores(data = counts_example, 
           metadata = metadata_example, 
           gene_sets = list(Senescence=SimpleSenescenceSignature_bidirectional),
           ColorVariable = "SenescentType", 
           GroupingVariable="Condition",  
           method ="ssGSEA", 
           ColorValues = senescence_triggers_colors, 
           ConnectGroups=TRUE, 
           ncol = NULL, 
           nrow = NULL, 
           widthTitle=24, 
           limits = NULL, 
           legend_nrow = 1, 
           pointSize=4,
           cond_cohend=cond_cohend,
           title="Marthandan et al. 2016",
           labsize=7, 
           titlesize = 10)  
```

<img src="man/figures/README-examplessGSEA_bi-1.png" width="40%" />

#### Ranking method

The following example uses the **“ranking”** method for score
calculation, both for unidirectional and bidirectional signatures.

``` r
 
# df_Scores <- CalculateScores(data = counts_example, 
#                              metadata = metadata_example, 
#                              method = "ranking", 
#                              gene_sets = list(Senescence=SimpleSenescenceSignature))

senescence_triggers_colors <- c(
  "none" = "#E57373",  # Soft red   
  "Telomere shortening" = "#4FC3F7"  # Vivid sky blue  
)

cond_cohend <- list(A=c("Senescent"),  
                    B=c("Proliferative"))

PlotScores(data = counts_example, 
           metadata = metadata_example, 
           gene_sets = list(Senescence=SimpleSenescenceSignature),
           ColorVariable = "SenescentType", 
           GroupingVariable="Condition",  
           method ="ranking", 
           ColorValues = senescence_triggers_colors, 
           ConnectGroups=TRUE, 
           ncol = NULL, 
           nrow = NULL, 
           widthTitle=24, 
           limits = NULL, 
           legend_nrow = 1, 
           pointSize=4,
           cond_cohend=cond_cohend,
           title="Marthandan et al. 2016",
           labsize=7, 
           titlesize = 10)  
```

<img src="man/figures/README-ranking_unidirect-1.png" width="40%" />

``` r
 
# df_Scores <- CalculateScores(data = counts_example, 
#                              metadata = metadata_example, 
#                              method = "ranking", 
#                              gene_sets = list(Senescence=SimpleSenescenceSignature_bidirectional))

senescence_triggers_colors <- c(
  "none" = "#E57373",  # Soft red   
  "Telomere shortening" = "#4FC3F7"  # Vivid sky blue  
)

cond_cohend <- list(A=c("Senescent"),  
                    B=c("Proliferative"))

PlotScores(data = counts_example, 
           metadata = metadata_example, 
           gene_sets = list(Senescence=SimpleSenescenceSignature_bidirectional),
           ColorVariable = "SenescentType", 
           GroupingVariable="Condition",  
           method ="ranking", 
           ColorValues = senescence_triggers_colors, 
           ConnectGroups=TRUE, 
           ncol = NULL, 
           nrow = NULL, 
           widthTitle=24, 
           limits = NULL, 
           legend_nrow = 1, 
           pointSize=4,
           cond_cohend=cond_cohend,
           title="Marthandan et al. 2016",
           labsize=7, 
           titlesize = 10)  
```

<img src="man/figures/README-ranking_bidirect-1.png" width="40%" />

#### All methods

To compare various metrics across different condition combinations,
violin plots may not always be the best choice. In such cases, users can
set to generate a summary heatmap. The function will return one heatmap
per gene set, with rows corresponding to all possible combinations of
values in the .

``` r
 
PlotScores(data = counts_example, 
           metadata = metadata_example,  
           gene_sets=list(Senescence_Bidirectional = SimpleSenescenceSignature_bidirectional,
                          Senescence  = SimpleSenescenceSignature), 
           GroupingVariable="Condition",  
           method ="all",   
           ncol = NULL, 
           nrow = NULL, 
           widthTitle=30, 
           limits = NULL,   
           title="Marthandan et al. 2016", 
           titlesize = 12,
           ColorValues = NULL)  
```

<img src="man/figures/README-heatmap_all-1.png" width="80%" />

### False Discovery Rate (FDR) Calculations

#### Simulation Based Methods

The user can assess the significance of gene signature scores by
comparing observed effect sizes against those originated by random
signatures. For each original gene signature, the function calculates
the observed Cohen’s d (and p‑value) using (`GroupingVariable`). It then
generates a number of simulated signatures (`number_of_sims`) by
randomly sampling genes from a user provided gene list (`gene_list`) and
computes their Cohen’s d values. The simulation results are visualised
as violin plots that display the distribution of Cohen’s d values for
each method, overlaid with the observed values of the original
signatures, and a 95th percentile threshold. Significance is indicated
by distinct point shapes based on the associated p‑value.

``` r
 
FDR_Simulation(data = counts_example,
               metadata = metadata_example,
               original_signatures = list(Senescence_Bidirectional = SimpleSenescenceSignature_bidirectional,
                          Senescence  = SimpleSenescenceSignature),
               gene_list = row.names(counts_example),
               number_of_sims = 100,
               title_for_plot = "Marthandan et al. 2016",
               GroupingVariable = "Condition")
```

<img src="man/figures/README-FDRSim-1.png" width="80%" />
