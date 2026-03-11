# Calculate Gene Signature Scores using Ranking Approach

Computes gene signature scores for each sample by ranking the expression
of signature genes in the dataset and normalizing the score based on the
total number of genes.

## Usage

``` r
CalculateScores_Ranking(data, metadata = NULL, gene_sets)
```

## Arguments

- data:

  A data frame where rows represent genes, columns represent samples,
  and values correspond to gene expression levels. **(Required)**

- metadata:

  A data frame containing sample metadata. The first column must contain
  sample names. **(Optional)**

- gene_sets:

  A named list of gene sets. **(Required)** For unidirectional gene
  sets, provide a named list where each element is a vector of gene
  names. For bidirectional gene sets, provide a named list where each
  element is a data frame with two columns:

  - The first column: gene names.

  - The second column: expected direction (1 for upregulated, -1 for
    downregulated).

## Value

A named list of data frames, where each data frame contains:

- `sample`: Sample name.

- `score`: Normalized ranking score for the given gene signature.

- Additional metadata columns (if `metadata` is provided).

## Details

- The function first validates inputs and extracts relevant genes from
  the dataset.

- For **unidirectional** signatures, it computes rankings based on gene
  expression levels.

- For **bidirectional** signatures, it computes separate rankings for
  upregulated and downregulated genes, then calculates a final score by
  subtracting downregulated rankings from upregulated rankings.

- The final scores are normalized by dividing by the total number of
  genes.

- This metric is not suitable to compare absolute values between
  different gene sets, i.e. should be used only for relative comparisons
  between samples when using the same gene set.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example dataset with 5 genes (rows) and 3 samples (columns)
set.seed(123)
data <- as.data.frame(matrix(runif(15, 1, 100), nrow = 5, ncol = 3))
rownames(data) <- paste0("Gene_", 1:5)
colnames(data) <- paste0("Sample_", 1:3)

# Unidirectional gene set example
gene_sets <- list(Signature1 = c("Gene_1", "Gene_3", "Gene_5"))

# Compute scores
scores <- CalculateScores_Ranking(data, gene_sets = gene_sets)
print(scores)
} # }
```
