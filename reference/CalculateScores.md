# Calculate Gene Signature Scores using Score-Based Approaches

This function calculates a gene signature score for each sample based on
one or more predefined gene sets (signatures).

## Usage

``` r
CalculateScores(
  data,
  metadata,
  gene_sets,
  method = c("ssGSEA", "logmedian", "ranking", "all")
)
```

## Arguments

- data:

  A data frame of normalized (non-transformed) counts where each row is
  a gene and each column is a sample. The row names should contain gene
  names, and the column names should contain sample identifiers.
  **(Required)**

- metadata:

  A data frame describing the attributes of each sample. Each row
  corresponds to a sample and each column to an attribute. The first
  column of `metadata` should be the sample identifiers (i.e., the
  column names of `data`). Defaults to `NULL` if no metadata is
  provided.

- gene_sets:

  Gene set input. **(Required)**

  If using **unidirectional** gene sets, provide a named list where each
  element is a vector of gene names representing a gene signature. The
  names of the list elements should correspond to the labels for each
  signature.

  If using **bidirectional** gene sets, provide a named list where each
  element is a data frame. The names of the list elements should
  correspond to the labels for each signature, and each data frame
  should contain the following structure:

  - The **first column** should contain gene names.

  - The **second column** should indicate the expected direction of
    enrichment (1 for upregulated genes, -1 for downregulated genes).

- method:

  A character string indicating the scoring method to use. Options are
  `"ssGSEA"`, `"logmedian"`, `"ranking"`, or `"all"` (to compute scores
  using all methods). Defaults to `"logmedian"`.

## Value

If a single method is chosen, a data frame containing the calculated
scores for each gene signature, including metadata if provided. If
`method = "all"`, a list is returned where each element corresponds to a
scoring method and contains the respective data frame of scores.

- sample:

  The sample identifier (matching the column names of the input data).

- score:

  The calculated gene signature score for the corresponding sample.

- (metadata):

  Any additional columns from the `metadata` data frame provided by the
  user, if available.

## Details

- `ssGSEA`:

  Uses the single-sample Gene Set Enrichment Analysis (ssGSEA) method to
  compute an enrichment score for each signature in each sample. This
  method uses an adaptation from the the `gsva()` function from the
  `GSVA` package to compute an enrichment score, representing the
  absolute enrichment of each gene set in each sample.

- `logmedian`:

  Computes, for each sample, the score as the sum of the normalized
  (log2-median-centered) expression values of the signature genes
  divided by the number of genes in the signature.

- `ranking`:

  Computes gene signature scores for each sample by ranking the
  expression of signature genes in the dataset and normalizing the score
  based on the total number of genes.

- `all`:

  Computes gene signature scores using all three methods (`ssGSEA`,
  `logmedian`, and `ranking`). The function returns a list containing
  the results of each method.

## Examples

``` r
# Simulate positive gene expression data (genes as rows, samples as columns)
set.seed(42)
expr <- as.data.frame(matrix(rexp(60, rate = 0.2), nrow = 6, ncol = 10))
rownames(expr) <- paste0("Gene", 1:6)
colnames(expr) <- paste0("Sample", 1:10)

# Simulate metadata for samples
metadata <- data.frame(
  sample = colnames(expr),
  Group = rep(c("A", "B"), each = 5)
)

# Define two simple gene sets
gene_sets <- list(
  Signature1 = c("Gene1", "Gene2", "Gene3"),
  Signature2 = c("Gene4", "Gene5", "Gene6")
)

# Calculate logmedian scores
scores_logmedian <- CalculateScores(
  data = expr,
  metadata = metadata,
  gene_sets = gene_sets,
  method = "logmedian"
)
#> Considering unidirectional gene signature mode for signature Signature1
#> Considering unidirectional gene signature mode for signature Signature2
head(scores_logmedian)
#> $Signature1
#>      sample       score Group
#> 1   Sample1 -0.59042384     A
#> 2  Sample10 -0.16112441     B
#> 3   Sample2 -0.12502535     A
#> 4   Sample3 -0.78459915     A
#> 5   Sample4  0.96089184     A
#> 6   Sample5  0.59656555     A
#> 7   Sample6  0.53661412     B
#> 8   Sample7 -0.97350624     B
#> 9   Sample8  0.50957254     B
#> 10  Sample9  0.01641426     B
#> 
#> $Signature2
#>      sample       score Group
#> 1   Sample1 -0.62074224     A
#> 2  Sample10  0.80397814     B
#> 3   Sample2  0.64215946     A
#> 4   Sample3 -0.59233078     A
#> 5   Sample4  0.32478963     A
#> 6   Sample5 -0.55206694     A
#> 7   Sample6  0.02289467     B
#> 8   Sample7  0.28685526     B
#> 9   Sample8  0.66827204     B
#> 10  Sample9  0.10551493     B
#> 

# Calculate all score types
scores_all <- CalculateScores(
  data = expr,
  metadata = metadata,
  gene_sets = gene_sets,
  method = "all"
)
#> Considering unidirectional gene signature mode for signature Signature1
#> No id variables; using all as measure variables
#> Considering unidirectional gene signature mode for signature Signature2
#> No id variables; using all as measure variables
#> Considering unidirectional gene signature mode for signature Signature1
#> Considering unidirectional gene signature mode for signature Signature2
#> Considering unidirectional gene signature mode for signature Signature1
#> Considering unidirectional gene signature mode for signature Signature2
lapply(scores_all, head)
#> $ssGSEA
#> $ssGSEA$Signature1
#>      sample        score Group
#> 1   Sample1 -0.036060828     A
#> 2  Sample10 -0.123345243     B
#> 3   Sample2 -0.359875397     A
#> 4   Sample3 -0.099807778     A
#> 5   Sample4  0.075203365     A
#> 6   Sample5  0.292611744     A
#> 7   Sample6  0.088004298     B
#> 8   Sample7 -0.484971533     B
#> 9   Sample8 -0.008588505     B
#> 10  Sample9 -0.099807778     B
#> 
#> $ssGSEA$Signature2
#>      sample        score Group
#> 1   Sample1  0.116563042     A
#> 2  Sample10  0.196723601     B
#> 3   Sample2  0.403808409     A
#> 4   Sample3  0.173737585     A
#> 5   Sample4  0.004629832     A
#> 6   Sample5 -0.231201575     A
#> 7   Sample6 -0.008588505     B
#> 8   Sample7  0.505619476     B
#> 9   Sample8  0.088004298     B
#> 10  Sample9  0.173737585     B
#> 
#> 
#> $logmedian
#> $logmedian$Signature1
#>      sample       score Group
#> 1   Sample1 -0.59042384     A
#> 2  Sample10 -0.16112441     B
#> 3   Sample2 -0.12502535     A
#> 4   Sample3 -0.78459915     A
#> 5   Sample4  0.96089184     A
#> 6   Sample5  0.59656555     A
#> 7   Sample6  0.53661412     B
#> 8   Sample7 -0.97350624     B
#> 9   Sample8  0.50957254     B
#> 10  Sample9  0.01641426     B
#> 
#> $logmedian$Signature2
#>      sample       score Group
#> 1   Sample1 -0.62074224     A
#> 2  Sample10  0.80397814     B
#> 3   Sample2  0.64215946     A
#> 4   Sample3 -0.59233078     A
#> 5   Sample4  0.32478963     A
#> 6   Sample5 -0.55206694     A
#> 7   Sample6  0.02289467     B
#> 8   Sample7  0.28685526     B
#> 9   Sample8  0.66827204     B
#> 10  Sample9  0.10551493     B
#> 
#> 
#> $ranking
#> $ranking$Signature1
#>      sample    score Group
#> 1   Sample1 1.666667     A
#> 2  Sample10 1.500000     B
#> 3   Sample2 1.166667     A
#> 4   Sample3 1.500000     A
#> 5   Sample4 1.833333     A
#> 6   Sample5 2.166667     A
#> 7   Sample6 1.833333     B
#> 8   Sample7 1.000000     B
#> 9   Sample8 1.666667     B
#> 10  Sample9 1.500000     B
#> 
#> $ranking$Signature2
#>      sample    score Group
#> 1   Sample1 1.833333     A
#> 2  Sample10 2.000000     B
#> 3   Sample2 2.333333     A
#> 4   Sample3 2.000000     A
#> 5   Sample4 1.666667     A
#> 6   Sample5 1.333333     A
#> 7   Sample6 1.666667     B
#> 8   Sample7 2.500000     B
#> 9   Sample8 1.833333     B
#> 10  Sample9 2.000000     B
#> 
#> 
```
