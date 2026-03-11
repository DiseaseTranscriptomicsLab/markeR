# Plot Gene Set Scores by Group or Continuous Variable

This function computes and visualizes gene set enrichment scores using
various methods, optionally comparing across groups or numeric
variables. It supports categorical and numeric comparisons, statistical
testing, Cohen's d effect sizes, and visualizations such as heatmaps and
volcano plots.

## Usage

``` r
PlotScores_Categorical(
  data,
  metadata,
  gene_sets,
  method = c("ssGSEA", "logmedian", "ranking"),
  ColorVariable = NULL,
  GroupingVariable = NULL,
  ColorValues = NULL,
  ConnectGroups = FALSE,
  ncol = NULL,
  nrow = NULL,
  title = NULL,
  widthTitle = 10,
  titlesize = 12,
  limits = NULL,
  legend_nrow = NULL,
  pointSize = 2,
  xlab = NULL,
  labsize = 10,
  compute_cohen = TRUE,
  cond_cohend = NULL,
  pvalcalc = FALSE,
  mode = c("simple", "medium", "extensive"),
  widthlegend = 22,
  cohen_threshold = 0.6,
  colorPalette = "Set3"
)
```

## Arguments

- data:

  A data frame of Normalised (non-transformed) counts where each row is
  a gene and each column is a sample. Row names should contain gene
  names, and column names should contain sample identifiers.
  **(Required)**

- metadata:

  A data frame describing the attributes of each sample, where each row
  corresponds to a sample and each column to an attribute. The first
  column should contain sample identifiers (i.e., the column names of
  `data`). **(Required if method = "all")**

- gene_sets:

  Gene set input. **(Required)**

  - **Unidirectional gene sets**: Provide a named list where each
    element is a vector of gene names representing a gene signature.

  - **Bidirectional gene sets**: Provide a named list where each element
    is a data frame with two columns:

    - The **first column** contains gene names.

    - The **second column** indicates the expected direction of
      enrichment (1 for upregulated genes, -1 for downregulated genes).

- method:

  A character string indicating the scoring method to use. Options are
  `"ssGSEA"`, `"logmedian"` or `"ranking"`. Defaults to `"logmedian"`.

- ColorVariable:

  Optional. Name of the metadata column to use for point color in plots.

- GroupingVariable:

  Optional. Name of the metadata column to use for group comparison.

- ColorValues:

  Optional. Named vector of colors to use for each group in
  `ColorVariable` or `GroupingVariable`.

- ConnectGroups:

  Logical. If TRUE, connects points of the same sample across
  conditions.

- ncol:

  Number of columns in the facet layout of the plot.

- nrow:

  Number of rows in the facet layout of the plot.

- title:

  Optional. Main title of the plot.

- widthTitle:

  Numeric. Width of the title area (for alignment purposes).

- titlesize:

  Numeric. Font size of the title text.

- limits:

  Optional numeric vector of length 2 specifying y-axis limits.

- legend_nrow:

  Optional. Number of rows in the plot legend.

- pointSize:

  Numeric. Size of the points in the plots.

- xlab:

  Optional. Label for the x-axis.

- labsize:

  Numeric. Font size for axis and facet labels.

- compute_cohen:

  Logical. If TRUE, computes Cohen's d effect sizes between groups.

- cond_cohend:

  Optional. Specify a condition or comparison subset for calculating
  Cohen's d.

- pvalcalc:

  Logical. If TRUE, computes p-values for group comparisons.

- mode:

  Character string indicating comparison complexity. Options:
  `"simple"`, `"medium"`, `"extensive"`.

- widthlegend:

  Numeric. Width of the legend area in volcano plots.

- cohen_threshold:

  Numeric. Cohen's d threshold to highlight effect size in volcano plots
  (default = 0.6).

- colorPalette:

  Character. Name of RColorBrewer palette for coloring (default =
  `"Set3"`).

## Value

A `ggplot` or a
[`ggpubr::ggarrange`](https://rpkgs.datanovia.com/ggpubr/reference/ggarrange.html)
object depending on the input and parameters:

- If `GroupingVariable` is `NULL`, returns a faceted grid of density
  plots (one per gene set).

- If `GroupingVariable` is provided and `method != "all"`, returns a
  faceted grid of violin plots overlaid with jittered sample points and
  median bars, optionally annotated with Cohen's d or f and p-values.

- Each individual plot corresponds to one gene set score computed using
  the selected method.

## Details

Four methods are available:

- **ssGSEA**: Uses the single-sample Gene Set Enrichment Analysis
  (ssGSEA) method to compute an enrichment score for each signature in
  each sample using an adaptation of the `gsva()` function from the
  `GSVA` package.

- **logmedian**: Computes the score as the sum of the Normalised
  (log2-median-centered) expression values of the signature genes
  divided by the number of genes in the signature.

- **ranking**: Computes gene signature scores for each sample by ranking
  the expression of signature genes in the dataset and normalizing the
  score based on the total number of genes.

- **all**: Computes gene signature scores using all three methods
  (`ssGSEA`, `logmedian`, and `ranking`). Returns a heatmap summarizing
  Cohen's d for all metric combinations of the variables of interest.

Depending on the method and the type of variable (categorical, numeric,
or `NULL`), the function produces different plots:

- **If `method = "all"`** and the variable is **categorical**, a heatmap
  of Cohen's d or F statistics and a volcano plot showing contrasts
  between all groups of that variable are produced.

- **If `method = "all"`** and the variable is **numeric**, a heatmap of
  Cohen's f and a volcano plot are produced.

- **If `method != "all"`** and the variable is **categorical**, a violin
  plot for each signature is generated.

- **If `method != "all"`** and the variable is `NULL`, a density plot of
  the score distribution is displayed.

- **If `method != "all"`** and the variable is **numeric**, a scatter
  plot is created to show the relationship between the scores and the
  numeric variable.

- **If `method = "all"`** and the variable is **categorical**, the
  function returns a heatmap of Cohen's d or F statistics and a volcano
  plot showing contrasts between all groups of that variable.

- **If `method = "all"`** and the variable is **numeric**, a heatmap of
  Cohen's f and a volcano plot will be produced.

- **If `method != "all"`** and the variable is **categorical**, a violin
  plot for each signature will be displayed.

- **If `method != "all"`** and the variable is `NULL`, a density plot of
  the score distribution will be displayed.

- **If `method != "all"`** and the variable is **numeric**, a scatter
  plot will be generated to show the relationship between the scores and
  the numeric variable.
