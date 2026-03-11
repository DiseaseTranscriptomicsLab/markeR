# Volcano Plot of Cohen's d Effect Sizes and Adjusted p-values

This function computes Cohen's d effect sizes and adjusted p-values for
multiple gene signatures across defined contrasts, and generates a
volcano plot (Cohen's d vs -log10(padj)) using `ggplot2`. Each point
represents a method-signature pair, faceted by contrast.

## Usage

``` r
Volcano_Cohen(
  cohenlist,
  titlesize = 12,
  ColorValues = NULL,
  title = NULL,
  widthlegend = 22,
  pointSize = 3,
  sig_threshold = 0.05,
  cohen_threshold = 0.5,
  colorPalette = "Set3",
  nrow = NULL,
  ncol = NULL
)
```

## Arguments

- cohenlist:

  A named list from `CohenD_allConditions`. Each element is a list with:

  - `CohenD`: A data frame where rows are methods and columns are group
    contrasts (formatted as "Group1:Group2"), containing the computed
    Cohen's d effect sizes.

  - `PValue`: A data frame with the same structure as `CohenD`
    containing the corresponding p-values.

  - `padj`: A data frame with the same structure as `PValue` containing
    the corresponding p-values corrected using the BH method, for all
    signatures and contrasts, and by method.

- titlesize:

  Integer. Size of the facet strip titles. Default is 12.

- ColorValues:

  Character vector of colors used to distinguish signatures. If NULL,
  colors are automatically generated.

- title:

  Optional title for the overall plot.

- widthlegend:

  Integer. Width used to wrap long signature names. Default is 22.

- pointSize:

  Numeric. Size of the points in the plot. Default is 3.

- sig_threshold:

  Numeric. Adjusted p-value threshold for significance. Default is 0.05.

- cohen_threshold:

  Numeric. Effect size threshold. Default is 0.5.

- colorPalette:

  Character. Name of RColorBrewer palette to use if `ColorValues` is not
  provided. Default is "Set3".

- nrow:

  Optional numeric value specifying the number of rows in the grid
  layout. If `NULL`, a near-square grid is computed.

- ncol:

  Optional numeric value specifying the number of columns in the ggplot
  facet. If `NULL`, a near-square grid is computed.

## Value

A `ggplot` object showing a faceted volcano plot of Cohen's d effect
sizes across signatures and methods for each contrast.

## See also

[`CohenD_allConditions`](https://diseasetranscriptomicslab.github.io/markeR/reference/CohenD_allConditions.md)
