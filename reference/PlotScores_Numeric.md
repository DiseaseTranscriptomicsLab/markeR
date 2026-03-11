# Plot Gene Signature Scores with Continuous Variables

This function visualizes gene signature scores using scatter plots and
regression lines across a continuous metadata variable. Signature scores
are computed per sample using one of three methods: `"ssGSEA"`,
`"logmedian"`, or `"ranking"`. Optionally, the effect size (Cohen's f)
and p-value for the association between the signature score and the
continuous variable can be computed and displayed.

## Usage

``` r
PlotScores_Numeric(
  data,
  metadata,
  gene_sets,
  method = c("ssGSEA", "logmedian", "ranking"),
  Variable = NULL,
  ColorValues = NULL,
  ncol = NULL,
  nrow = NULL,
  title = NULL,
  widthTitle = 10,
  titlesize = 12,
  limits = NULL,
  pointSize = 2,
  xlab = NULL,
  labsize = 10,
  compute_cohen = TRUE,
  pvalcalc = FALSE,
  colorPalette = "Set3",
  cor = c("pearson", "spearman", "kendall")
)
```

## Arguments

- data:

  A data frame of Normalised (non-transformed) gene expression counts.
  Rows are genes, columns are samples. Row names should be gene names,
  and column names should match sample identifiers in `metadata`.

- metadata:

  A data frame where each row corresponds to a sample and contains
  sample-level attributes (e.g., clinical or experimental metadata).
  Must include a column matching the sample IDs in `data`.

- gene_sets:

  A list of gene sets (signatures). Each element is either a character
  vector of gene names or a data frame with gene names and enrichment
  direction (1 for upregulated, -1 for downregulated).

- method:

  Scoring method to use. One of `"ssGSEA"`, `"logmedian"`, or
  `"ranking"`. Default is `"logmedian"`.

- Variable:

  Name of the continuous variable in `metadata` to use on the x-axis for
  scoring association.

- ColorValues:

  (Optional) A named vector defining the color for the plotted points.
  If NULL, defaults to a preset color.

- ncol, nrow:

  Number of columns and rows in the facet grid layout. If NULL, computed
  automatically.

- title:

  Optional string for the overall title of the plot grid.

- widthTitle:

  Maximum character width for titles before inserting line breaks.
  Default is 10.

- titlesize:

  Numeric value for the font size of plot titles. Default is 12.

- limits:

  Optional numeric vector of length 2 to define y-axis limits.

- pointSize:

  Size of the plotted points. Default is 2.

- xlab:

  Optional label for the x-axis. If NULL, defaults to the name of
  `Variable`.

- labsize:

  Font size for axis labels. Default is 10.

- compute_cohen:

  Logical. If TRUE (default), computes Cohen's f effect size for the
  association between signature score and the continuous variable.

- pvalcalc:

  Logical. If TRUE, includes the p-value in the plot subtitle. Default
  is FALSE.

- colorPalette:

  Name of the RColorBrewer palette for coloring. Default is "Set3".
  Currently unused but kept for consistency.

- cor:

  Character string indicating the correlation method to be used in
  [`ggpubr::stat_cor()`](https://rpkgs.datanovia.com/ggpubr/reference/stat_cor.html).
  Options are "pearson" (default), "kendall", or "spearman".

## Value

A ggplot2 object or a multi-plot figure showing scatter plots for each
gene signature, with linear regression lines and optional statistical
annotations.

## Details

For each gene signature, the function:

- Computes a signature score per sample using the selected method.

- Plots the score against a continuous metadata variable (`Variable`).

- Adds a regression line and optionally computes and displays Cohen's f
  effect size and p-value.

- Returns a faceted grid of ggplots, arranged by `ncol` and `nrow`.

This version of the function is specifically tailored for use with
continuous variables.
