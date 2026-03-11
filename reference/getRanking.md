# Get Gene Expression Ranking

Computes the rank sum of a given gene set within a sample based on its
expression level.

## Usage

``` r
getRanking(data, sample, geneset)
```

## Arguments

- data:

  A data frame where rows represent genes, columns represent samples,
  and values correspond to expression levels.

- sample:

  A character string specifying the sample name (column in `data`).

- geneset:

  A vector of gene names to be ranked.

## Value

The sum of the ranks of the genes found in the sample.

## Details

- The function orders gene expression levels from lowest to highest.

- It then determines the rank of each gene in `geneset` and returns the
  sum of these ranks.

- If some genes are missing, they are omitted from the ranking
  calculation.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example dataset with 5 genes and 3 samples
set.seed(123)
data <- as.data.frame(matrix(runif(15, 1, 100), nrow = 5, ncol = 3))
rownames(data) <- paste0("Gene_", 1:5)
colnames(data) <- paste0("Sample_", 1:3)

# Define gene set
geneset <- c("Gene_1", "Gene_3", "Gene_5")

# Compute ranking for Sample_1
rank_score <- getRanking(data, "Sample_1", geneset)
print(rank_score)
} # }
```
