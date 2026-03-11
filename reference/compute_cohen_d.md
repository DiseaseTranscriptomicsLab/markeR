# Compute Pairwise Cohen\\s d and P-Values

Computes Cohen\\s d effect sizes and corresponding p-values for all
pairwise comparisons of a grouping variable in a data frame.

## Usage

``` r
compute_cohen_d(
  dfScore,
  variable,
  quantitative_var = "score",
  mode = c("simple", "medium", "extensive")
)
```

## Arguments

- dfScore:

  A data frame containing at least one numeric column and a grouping
  variable. Output from flatten_results.

- variable:

  A string specifying the name of the categorical grouping column in
  `dfScore`.

- quantitative_var:

  A string specifying the name of the numeric column (default is
  `"score"`).

## Value

A data frame with the following columns:

- Group1:

  The first group in the pair.

- Group2:

  The second group in the pair.

- CohenD:

  The computed Cohen\\s d effect size for the comparison.

- PValue:

  The p-value from a t-test comparing the two groups.
