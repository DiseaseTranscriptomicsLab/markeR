# Calculate ssGSEA Scores for Bidirectional Gene Signatures

Computes single-sample Gene Set Enrichment Analysis (ssGSEA) scores for
each sample using a bidirectional gene signature (separating upregulated
and downregulated genes).

## Usage

``` r
CalculateScores_ssGSEA_bidirectional(data, signature)
```

## Arguments

- data:

  A data frame of normalized (non-transformed) counts where rows are
  genes and columns are samples.

- signature:

  A data frame with:

  - The **first column** containing gene names.

  - The **second column** (`Signal`) indicating the expected direction
    of enrichment (1 for upregulated genes, -1 for downregulated genes).

## Value

A data frame containing:

- `sample`: Sample name.

- `score`: Final ssGSEA enrichment score (computed as the difference
  between upregulated and downregulated scores).

## Details

- The input gene expression matrix (`data`) is log2-transformed before
  applying ssGSEA.

- Upregulated and downregulated genes are analyzed separately.

- As both upregulated and downregulated genes are present, the final
  score is computed as: \$\$score = (score\_{up} \\
  \frac{\|up\\genes\|}{\|total\\genes\|}) - (score\_{down} \\
  \frac{\|down\\genes\|}{\|total\\genes\|})\$\$

- If no downregulated genes are present, only the upregulated score is
  used.

- The results are reshaped into a long-format data frame with one score
  per sample.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example dataset with 5 genes (rows) and 3 samples (columns)
set.seed(123)
data <- matrix(runif(15, 1, 100), nrow = 5, ncol = 3)
rownames(data) <- paste0("Gene_", 1:5)
colnames(data) <- paste0("Sample_", 1:3)

# Define a bidirectional gene signature
signature <- data.frame(Gene = c("Gene_1", "Gene_3", "Gene_5"),
                        Signal = c(1, -1, 1))

# Compute scores
scores <- CalculateScores_ssGSEA_bidirectional(data, signature = signature)
print(scores)
} # }
```
