# GSEA Variable Association

This function assesses the association between gene expression (or
another molecular score) and metadata variables using differential
expression (DE) analysis and Gene Set Enrichment Analysis (GSEA). It
generates all possible contrasts for categorical variables and uses
linear modeling for continuous variables.

## Usage

``` r
GSEA_VariableAssociation(
  data,
  metadata,
  cols,
  stat = NULL,
  mode = c("simple", "medium", "extensive"),
  gene_set,
  nonsignif_color = "grey",
  signif_color = "red",
  saturation_value = NULL,
  sig_threshold = 0.05,
  widthlabels = 18,
  labsize = 10,
  titlesize = 14,
  pointSize = 5,
  ignore_NAs = FALSE,
  printplt = TRUE,
  p.adjust.method = "BH"
)
```

## Arguments

- data:

  A matrix or data frame containing gene expression data, where rows
  represent genes and columns represent samples.

- metadata:

  A data frame containing sample metadata with at least one column
  corresponding to the variables of interest.

- cols:

  A character vector specifying the metadata columns (variables) to
  analyse.

- stat:

  Optional. The statistic to use for ranking genes before GSEA. If
  `NULL`, it is automatically determined based on the gene set:

  - `"B"` for gene sets with **no known direction** (vectors).

  - `"t"` for **unidirectional** or **bidirectional** gene sets (data
    frames).

  - If provided, this argument overrides the automatic selection.

- mode:

  A string specifying the level of detail for contrasts. Options are:

  - `"simple"`: Performs the minimal number of pairwise comparisons
    between individual group levels (e.g., A - B, A - C). Default.

  - `"medium"`: Includes comparisons between one group and the union of
    all other groups (e.g., A - (B + C + D)), enabling broader contrasts
    beyond simple pairs.

  - `"extensive"`: Allows for all possible algebraic combinations of
    group levels (e.g., (A + B) - (C + D)), supporting flexible and
    complex contrast definitions.

- gene_set:

  A named list defining the gene sets for GSEA. **(Required)**

  - If using **unidirectional** gene sets, provide a list where each
    element is a vector of gene names representing a signature.

  - If using **bidirectional** gene sets, provide a list where each
    element is a data frame:

  - The **first column** should contain gene names.

  - The **second column** should indicate the expected direction of
    enrichment (`1` for upregulated, `-1` for downregulated).

- nonsignif_color:

  A string specifying the color for the middle of the adjusted p-value
  gradient. Default is `"white"`. Lower limit correspond to the value of
  `sig_threshold`.

- signif_color:

  A string specifying the color for the low end of the adjusted p-value
  gradient until the value chosen for significance (`sig_threshold`).
  Default is `"red"`.

- saturation_value:

  A numeric value specifying the lower limit of the adjusted p-value
  gradient, below which the color will correspond to `signif_color`.
  Default is the results' minimum, unless that value is above the
  sig_threshold; in that case, it is 0.001.

- sig_threshold:

  A numeric value specifying the threshold for significance
  visualization in the plot. Default: `0.05`.

- widthlabels:

  An integer controlling the maximum width of contrast labels before
  text wrapping. Default: `18`.

- labsize:

  An integer controlling the axis text size in the plot. Default: `10`.

- titlesize:

  An integer specifying the plot title size. Default: `14`.

- pointSize:

  Numeric. The size of points in the lollipop plot (default is 5).

- ignore_NAs:

  Boolean (default: FALSE). Whether to ignore NAs in the metadata when
  fitting the linear model. If TRUE, rows with any NAs will be removed
  before analysis, leading to a loss of data to be fitted in the model.

- printplt:

  Boolean specifying if plot is to be printed. Default: `TRUE`.

- p.adjust.method:

  Character string specifying the method to use for multiple testing
  correction. Must be one of `"BH"` (Benjamini-Hochberg, default),
  `"holm"`, `"hommel"`, `"bonferroni"`, `"BY"` (Benjamini-Yekutieli),
  `"fdr"`, or `"none"`. Passed to
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html).

## Value

A list with two elements:

- `data`: A data frame containing the GSEA results, including normalized
  enrichment scores (NES), adjusted p-values, and contrasts.

- `plot`: A ggplot2 object visualizing the GSEA results as a lollipop
  plot.
