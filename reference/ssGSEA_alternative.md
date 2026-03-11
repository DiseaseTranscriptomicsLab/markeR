# Alternative Implementation of Single-Sample Gene Set Enrichment Analysis (ssGSEA)

This function computes an enrichment score for each sample using an
alternative single-sample Gene Set Enrichment Analysis (ssGSEA) method.
It first maps gene sets to the gene indices present in the expression
matrix, then ranks the genes for each sample, and finally calculates a
weighted enrichment score based on the cumulative differences between
in-set and out-of-set gene ranks. Source:
https://rpubs.com/pranali018/SSGSEA

## Usage

``` r
ssGSEA_alternative(
  X,
  gene_sets,
  alpha = 0.25,
  scale = TRUE,
  norm = FALSE,
  single = TRUE
)
```

## Arguments

- X:

  A numeric matrix of gene expression values with rows representing
  genes and columns representing samples. Row names should correspond to
  gene identifiers.

- gene_sets:

  A list of gene sets, where each element is a vector of gene
  identifiers. The function will match these identifiers with the row
  names of `X`.

- alpha:

  A numeric value specifying the exponent used to weight the ranking
  scores. Default is `0.25`.

- scale:

  Logical; if `TRUE`, the cumulative difference is normalized by the
  total number of genes. Default is `TRUE`.

- norm:

  Logical; if `TRUE`, the enrichment scores are further normalized by
  the absolute difference between the maximum and minimum scores.
  Default is `FALSE`.

- single:

  Logical; if `TRUE`, the function returns the sum of the cumulative
  differences as the enrichment score. If `FALSE`, the maximum absolute
  cumulative difference is used. Default is `TRUE`.

## Value

A matrix of enrichment scores with rows corresponding to gene sets and
columns corresponding to samples.

## Details

The function performs the following steps:

1.  Maps each gene set to the indices of genes in `X` by matching gene
    identifiers.

2.  Computes column-wise rankings for the gene expression matrix using a
    ranking method (via the `colRanking` function) with tie resolution
    set to `'average'`.

3.  For each sample, orders the gene ranks in decreasing order.

4.  For each gene set in the sample, calculates:

    - The weighted contribution (`rank_alpha`) for genes in the set
      raised to the power of `alpha`.

    - The cumulative distribution functions (CDFs) for genes within the
      gene set (`step_cdf_pos`) and those not in the gene set
      (`step_cdf_neg`).

    - The difference between these CDFs, optionally scaled by the number
      of genes if `scale = TRUE`.

    - Depending on the `single` parameter, either the sum of the
      differences (if `TRUE`) or the maximum absolute difference (if
      `FALSE`) is used as the enrichment score for that gene set.

5.  Optionally normalizes the final enrichment scores by the range of
    values if `norm = TRUE`.

## Examples

``` r
if (FALSE) { # \dontrun{
  # Create a sample gene expression matrix:
  X <- matrix(rnorm(1000), nrow = 100, ncol = 10)
  rownames(X) <- paste0("gene", 1:100)

  # Define example gene sets:
  gene_sets <- list(
    set1 = sample(rownames(X), 10),
    set2 = sample(rownames(X), 15)
  )

  # Compute the ssGSEA enrichment scores:
  es <- ssGSEA_alternative(X, gene_sets, alpha = 0.25, scale = TRUE,
  norm = FALSE, single = TRUE)
  print(es)
} # }
```
