# Calculate Log-Median Scores for Unidirectional Gene Sets

Computes log-median-centered scores for gene signatures where all genes
are expected to be enriched in the same direction, or when direction is
not known.

## Usage

``` r
calculateScore_logmedian_unidirectional(data, signature)
```

## Arguments

- data:

  A data frame of normalized counts (genes as rows, samples as columns).

- signature:

  A vector of gene names or a data frame where the first column contains
  gene names.

## Value

A named vector with log-median-centered scores per sample.
