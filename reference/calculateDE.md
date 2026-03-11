# Calculate Differential Gene Expression Statistics using limma

This function computes differential gene expression statistics for each
gene using a linear model via the limma package. Users may supply a
custom design matrix directly via the `design` argument, or specify a
model formula (`lmexpression`) (e.g., `~0 + X` or `~X`) or variables
from `metadata` to build the design matrix. When contrasts are supplied,
they are applied using
[`limma::makeContrasts`](https://rdrr.io/pkg/limma/man/makeContrasts.html)
and
[`limma::contrasts.fit`](https://rdrr.io/pkg/limma/man/contrasts.fit.html).
Alternatively, when using `lmexpression` or a supplied `design`,
specific coefficient indices may be provided via `coefs` to extract the
corresponding gene-level statistics.

## Usage

``` r
calculateDE(
  data,
  metadata = NULL,
  variables = NULL,
  modelmat = NULL,
  contrasts = NULL,
  ignore_NAs = FALSE
)
```

## Arguments

- data:

  A numeric matrix of gene expression values with genes as rows and
  samples as columns. Row names must correspond to gene identifiers.
  Data should *not* be transformed (i.e., not log2 transformed).

- metadata:

  A data frame containing sample metadata used to build the design
  matrix (unless a design is provided directly).

- variables:

  A character vector specifying the variable(s) from `metadata` to use
  in the default linear model. Ignored if `lmexpression` or `design` is
  provided.

- modelmat:

  (Optional) A user-supplied design matrix. If provided, this design is
  used directly and `lmexpression` and `variables` are ignored. The
  order of samples in the design matrix should match the order in data.

- contrasts:

  A character vector specifying contrasts to be applied (e.g.,
  `c("A-B")`). If multiple contrasts are provided, the function returns
  a list of DE results (one per contrast). *Required* if `lmexpression`
  is NULL, optional otherwise. If not provided, the average expression
  profile of each condition will be returned instead of differential
  gene expression.

- ignore_NAs:

  Boolean (default: FALSE). Whether to ignore NAs in the metadata. If
  TRUE, rows with any NAs will be removed before analysis, leading to a
  loss of data to be fitted in the model. Only applicable if `variables`
  is provided.

## Value

A list of data-frames of differential expression statistics

## Details

The function fits a linear model with
[`limma::lmFit`](https://rdrr.io/pkg/limma/man/lmFit.html) and applies
empirical Bayes moderation with
[`limma::eBayes`](https://rdrr.io/pkg/limma/man/ebayes.html). Depending
on the input:

- If a design matrix is provided via `design`, that design is used
  directly.

- Otherwise, a design matrix is constructed using the `variables`
  argument (with no intercept).

- If contrasts are provided, they are applied using
  [`limma::makeContrasts`](https://rdrr.io/pkg/limma/man/makeContrasts.html)
  and
  [`limma::contrasts.fit`](https://rdrr.io/pkg/limma/man/contrasts.fit.html).

- If no contrasts are provided, the function returns all possible
  coefficients fitted in the linear model.

## Examples

``` r
# Simulate non-negative gene expression data (counts)
set.seed(123)
expr <- matrix(rpois(1000, lambda = 20), nrow = 100, ncol = 10)
rownames(expr) <- paste0("gene", 1:100)
colnames(expr) <- paste0("sample", 1:10)

# Simulate metadata with a group variable
metadata <- data.frame(
 sample = colnames(expr),
 Group = rep(c("A", "B"), each = 5)
)

# Differential expression for Group A vs Group B using variables
de_var <- calculateDE(
  data = expr,
  metadata = metadata,
  variables = "Group",
  contrasts = "A-B"
)
#> Using metadata column 'sample' to match samples (data column names).
head(de_var[["A-B"]])
#>              logFC  AveExpr         t      P.Value  adj.P.Val          B
#> gene99   0.8089957 4.271457  3.872148 0.0002292632 0.02292632  0.3805475
#> gene100  0.5616097 4.169432  2.789276 0.0066981874 0.33490937 -2.6692684
#> gene97   0.5475480 4.203882  2.573807 0.0120392169 0.33719450 -3.1834675
#> gene6    0.5345671 4.265823  2.530691 0.0134877800 0.33719450 -3.2822572
#> gene85   0.4461708 4.300741  2.194267 0.0313245769 0.57725317 -4.0041349
#> gene48  -0.4401609 4.364642 -2.090269 0.0399956221 0.57725317 -4.2091717

# Build equivalent design matrix manually
design <- model.matrix(~0 + Group, data = metadata)
colnames(design) <- c("A","B")

# Differential expression using the design matrix directly
de_mat <- calculateDE(
  data = expr,
  metadata = metadata,
  modelmat = design,
  contrasts = "A-B"
)
#> Using metadata column 'sample' to match samples (data column names).
head(de_mat[["A-B"]])
#>              logFC  AveExpr         t      P.Value  adj.P.Val          B
#> gene99   0.8089957 4.271457  3.872148 0.0002292632 0.02292632  0.3805475
#> gene100  0.5616097 4.169432  2.789276 0.0066981874 0.33490937 -2.6692684
#> gene97   0.5475480 4.203882  2.573807 0.0120392169 0.33719450 -3.1834675
#> gene6    0.5345671 4.265823  2.530691 0.0134877800 0.33719450 -3.2822572
#> gene85   0.4461708 4.300741  2.194267 0.0313245769 0.57725317 -4.0041349
#> gene48  -0.4401609 4.364642 -2.090269 0.0399956221 0.57725317 -4.2091717
```
