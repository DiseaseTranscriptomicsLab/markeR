# Compute ROC Curves and AUC Values for Gene Signature Scores

This function calculates Receiver Operating Characteristic (ROC) curves
and Area Under the Curve (AUC) values for gene signature scores across
different contrasts of a given categorical variable.

## Usage

``` r
ROCAUC_Scores_Calculate(
  data,
  metadata,
  gene_sets,
  method = c("logmedian", "ssGSEA", "ranking", "all"),
  variable,
  mode = c("simple", "medium", "extensive")
)
```

## Arguments

- data:

  A matrix or data frame of gene expression data (genes as rows, samples
  as columns).

- metadata:

  A data frame containing sample metadata, including the grouping
  variable.

- gene_sets:

  A named list of gene sets, where each entry is a character vector of
  gene names.

- method:

  A character string specifying the score calculation method. Options:
  `"logmedian"`, `"ssGSEA"`, `"ranking"`, or `"all"`.

- variable:

  A character string specifying the categorical variable for group
  comparisons.#'

- mode:

  A string specifying the level of detail for contrasts. Options are:

  - `"simple"`: Pairwise comparisons (e.g., A - B).

  - `"medium"`: Pairwise comparisons plus comparisons against the mean
    of other groups.

  - `"extensive"`: All possible groupwise contrasts, ensuring balance in
    the number of terms on each side.

## Value

A named list containing ROC curve data and AUC values for each method,
signature, and contrast.
