# Calculate Gene Signature Scores using ssGSEA

Computes an enrichment score for each gene signature in each sample
using the single-sample Gene Set Enrichment Analysis (ssGSEA).

## Usage

``` r
CalculateScores_ssGSEA(data, metadata = NULL, gene_sets)
```

## Arguments

- data:

  A data frame of normalized (non-transformed) counts where each row is
  a gene and each column is a sample.

- metadata:

  A data frame containing sample metadata (optional).

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

## Value

A list of data frames containing ssGSEA scores for each signature.
