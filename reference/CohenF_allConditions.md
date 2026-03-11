# Compute Cohen's f for All Gene Signatures Across a Categorical Variable

Computes Cohen's f effect sizes and corresponding p-values for all gene
signatures using scores calculated by multiple methods. The function
first computes gene signature scores using `CalculateScores` with the
"all" option, flattens the results, and fits linear models using the
specified variable to estimate effect sizes.

## Usage

``` r
CohenF_allConditions(
  data,
  metadata,
  gene_sets,
  variable,
  p.adjust.method = "BH"
)
```

## Arguments

- data:

  A data frame of gene expression data, with genes as rows and samples
  as columns.

- metadata:

  A data frame containing sample metadata. The first column should
  contain sample identifiers matching the column names of `data`.

- gene_sets:

  A named list of gene sets. For unidirectional gene sets, each element
  is a vector of gene names; for bidirectional gene sets, each element
  is a data frame where the first column contains gene names and the
  second column indicates the expected direction (1 for upregulated, -1
  for downregulated).

- variable:

  A string specifying the categorical variable in `metadata` used to
  model the gene signature scores.

- p.adjust.method:

  Character string specifying the method to use for multiple testing
  correction. Must be one of `"BH"` (Benjamini-Hochberg, default),
  `"holm"`, `"hommel"`, `"bonferroni"`, `"BY"` (Benjamini-Yekutieli),
  `"fdr"`, or `"none"`. Passed to
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html).

## Value

A named list where each element corresponds to a gene signature. Each
signature element is a list with three components:

- CohenF:

  A data frame where rows are scoring methods and columns are the
  variable used in the linear model (usually one column), containing the
  computed Cohen's f effect size.

- PValue:

  A data frame of the corresponding raw p-values from the linear model
  for each method.

- padj:

  A data frame of adjusted p-values (Benjamini-Hochberg method) across
  signatures and contrasts, per method.

## Details

This function is designed for use with categorical variables, where the
goal is to evaluate the overall group effect (e.g., using ANOVA) across
multiple levels.
