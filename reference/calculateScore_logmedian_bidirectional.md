# Calculate Log-Median Scores for Bidirectional Gene Sets

Computes gene signature scores considering both upregulated and
downregulated genes separately, then calculates a differential score by
subtracting downregulated from upregulated scores.

## Usage

``` r
calculateScore_logmedian_bidirectional(data, signature)
```

## Arguments

- data:

  A data frame of normalized counts (genes as rows, samples as columns).

- signature:

  A data frame with:

  - The **first column** containing gene names.

  - The **second column** specifying enrichment direction (`1` for
    upregulated, `-1` for downregulated).

## Value

A named vector with log-median-centered scores per sample.
