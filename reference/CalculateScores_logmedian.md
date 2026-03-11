# Calculate Gene Signature Scores using Log-Median Approach

Computes log2-median-centered scores for each sample based on gene
signature expression.

## Usage

``` r
CalculateScores_logmedian(data, metadata = NULL, gene_sets)
```

## Arguments

- data:

  A data frame of normalized counts where each row is a gene and each
  column is a sample.

- metadata:

  A data frame containing sample metadata (optional). If provided, the
  resulting scores will be merged with metadata.

- gene_sets:

  A named list representing gene sets. **(Required)**

  - **Unidirectional gene sets:** Each element should be a vector of
    gene names representing a signature. The names of the list elements
    serve as labels for the signatures.

  - **Bidirectional gene sets:** Each element should be a data frame
    with two columns:

    - First column: gene names.

    - Second column: expected direction of enrichment (1 for
      upregulated, -1 for downregulated).

## Value

A list of data frames containing log-median scores for each signature.
If `metadata` is provided, it is merged with the scores.

## Examples

``` r
if (FALSE) { # \dontrun{
data <- matrix(rnorm(1000), nrow = 100, ncol = 10)
colnames(data) <- paste0("Sample_", 1:10)
rownames(data) <- paste0("Gene_", 1:100)
gene_sets <- list(
  Signature_A = sample(rownames(data), 10),
  Signature_B = data.frame(Gene = sample(rownames(data), 10), Direction =
  sample(c(1, -1), 10, replace = TRUE))
)
scores <- CalculateScores_logmedian(data, gene_sets = gene_sets)
} # }
```
