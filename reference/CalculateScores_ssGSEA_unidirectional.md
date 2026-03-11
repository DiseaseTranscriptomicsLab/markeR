# Calculate ssGSEA Scores for Unidirectional Gene Signatures

Computes single-sample Gene Set Enrichment Analysis (ssGSEA) scores for
each sample using a unidirectional gene signature.

## Usage

``` r
CalculateScores_ssGSEA_unidirectional(data, signature)
```

## Arguments

- data:

  A data frame of normalized (non-transformed) counts where rows are
  genes and columns are samples.

- signature:

  A vector of gene names representing a unidirectional gene signature.

## Value

A data frame containing:

- `sample`: Sample name.

- `score`: ssGSEA enrichment score for the gene signature.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example dataset with 5 genes (rows) and 3 samples (columns)
set.seed(123)
data <- matrix(runif(15, 1, 100), nrow = 5, ncol = 3)
rownames(data) <- paste0("Gene_", 1:5)
colnames(data) <- paste0("Sample_", 1:3)

# Define a unidirectional gene signature
signature <- c("Gene_1", "Gene_3", "Gene_5")

# Compute scores
scores <- CalculateScores_ssGSEA_unidirectional(data, signature = signature)
print(scores)
} # }
```
