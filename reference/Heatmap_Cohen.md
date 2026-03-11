# Generate Heatmaps for Cohen's d Effect Sizes using ggplot2

This function computes Cohen's d effect sizes and corresponding p-values
for multiple gene signatures and produces individual heatmaps. Each
heatmap displays cell text showing the Cohen's d value along with its
p-value. The heatmaps are then arranged in a grid layout.

## Usage

``` r
Heatmap_Cohen(
  cohenlist,
  nrow = NULL,
  ncol = NULL,
  limits = NULL,
  widthTitle = 22,
  titlesize = 12,
  ColorValues = NULL,
  title = NULL
)
```

## Arguments

- cohenlist:

  A named list where each element corresponds to a gene signature.
  Output of `CohenD_allConditions`. Each signature element is a list
  with three components:

  CohenD

  :   A data frame where rows are methods and columns are group
      contrasts (formatted as "Group1:Group2"), containing the computed
      Cohen\\s d effect sizes.

  PValue

  :   A data frame with the same structure as `CohenD` containing the
      corresponding p-values.

  padj

  :   A data frame with the same structure as `PValue` containing the
      corresponding p-values corrected using the BH method, for all
      signatures and contrasts, and by method.

- nrow:

  Optional. An integer specifying the number of rows in the heatmap
  grid. If `NULL`, the number of rows is computed automatically.

- ncol:

  Optional. An integer specifying the number of columns in the heatmap
  grid. If `NULL`, the number of columns is computed automatically.

- limits:

  Optional. A numeric vector of length 2 specifying the color scale
  limits (e.g., `c(min, max)`). If `NULL`, the limits are determined
  from the data.

- widthTitle:

  An integer specifying the width used for wrapping gene set signature
  names in the heatmap titles. Default is 22.

- titlesize:

  An integer specifying the text size for each of the heatmap titles.
  Default is 12.

- ColorValues:

  A character vector specifying the colors for the gradient fill in the
  heatmaps. Default is `c("#F9F4AE", "#B44141")`.

- title:

  Title for the grid of plots.

## Value

A list with two elements:

- plt:

  A combined heatmap arranged in a grid using
  [`ggpubr::ggarrange`](https://rpkgs.datanovia.com/ggpubr/reference/ggarrange.html).

- data:

  A list containing the Cohen\\s d effect sizes and p-values for each
  gene signature, as computed by `CohenD_allConditions`.

## Details

The function first calculates Cohen\\s d effect sizes and corresponding
p-values for each gene signature using `CohenD_allConditions` (assumed
to be defined elsewhere in the package). The resulting matrices are
converted to a long format so that each cell in the heatmap can display
the Cohen\\s d value and its associated p-value (formatted as
`Cohen\'s d (p-value)`).

The heatmaps are then adjusted to display axis text and ticks only for
the left-most column and bottom row, and combined into a grid layout. If
neither `nrow` nor `ncol` are specified, the layout is automatically
determined to best approximate a square grid.

## See also

[`CohenD_allConditions`](https://diseasetranscriptomicslab.github.io/markeR/reference/CohenD_allConditions.md),
[`CohenF_allConditions`](https://diseasetranscriptomicslab.github.io/markeR/reference/CohenF_allConditions.md)
