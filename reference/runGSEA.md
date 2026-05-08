# Run Gene Set Enrichment Analysis (GSEA) for Multiple Contrasts

This function performs GSEA using `fgsea` for each contrast in a list of
differential expression results. It automatically determines the
appropriate ranking statistic based on the gene set format unless
specified by the user.

## Usage

``` r
runGSEA(
  DEGList,
  gene_sets,
  stat = NULL,
  ContrastCorrection = FALSE,
  nPermSimple = 10000,
  p.adjust.method = "BH"
)
```

## Arguments

- DEGList:

  A named list where each element represents a contrast and contains a
  data frame of differential expression results.

  - Each data frame must include at least the `"t"` statistic and the
    `"B"` statistic for each gene.

  - Row names should correspond to gene identifiers.

- gene_sets:

  A named list where each element represents a gene set. Each gene set
  can be:

  - A **vector of gene names** (for unidirectional gene sets).

  - A **data frame** with two columns:

    - Column 1: Gene names.

    - Column 2: Expected direction (`1` for upregulated genes, `-1` for
      downregulated genes).

- stat:

  Optional. The statistic to use for ranking genes before GSEA. If
  `NULL`, it is automatically determined based on the gene set:

  - `"B"` for gene sets with **no known direction** (vectors).

  - `"t"` for **unidirectional** or **bidirectional** gene sets (data
    frames).

  - If provided, this argument overrides the automatic selection.

- ContrastCorrection:

  Logical, default is `FALSE`. If `TRUE`, applies an additional multiple
  testing correction (Benjamini-Hochberg) across all contrasts returned
  in the `DEGList` results list. This accounts for the number of
  contrasts tested per signature and provides more stringent control of
  false discovery rate across multiple comparisons. If `FALSE`, the
  function only corrects for the number of gene sets.

- nPermSimple:

  Number of permutations in the simple fgsea implementation for
  preliminary estimation of P-values. Parameter from fgsea.

- p.adjust.method:

  Character string specifying the method to use for multiple testing
  correction. Must be one of `"BH"` (Benjamini-Hochberg, default),
  `"holm"`, `"hommel"`, `"bonferroni"`, `"BY"` (Benjamini-Yekutieli),
  `"fdr"`, or `"none"`. Passed to
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html).

## Value

A named list where each element corresponds to a contrast. Each contrast
contains a **single data frame** with GSEA results for all gene sets.
P-values are corrected for multiple testing based on all contrasts. The
result includes the standard `fgsea` output plus two additional columns:

- `pathway`: The name of the gene set.

- `stat_used`: The statistic used for ranking genes in that analysis
  (`"t"` or `"B"`).

## Examples

``` r
# Example input data
DEGList <- list(
  Contrast1 = data.frame(t = rnorm(100), B = rnorm(100), row.names = paste0("Gene", 1:100)),
  Contrast2 = data.frame(t = rnorm(100), B = rnorm(100), row.names = paste0("Gene", 1:100))
)

gene_sets <- list(
  UnidirectionalSet = c("Gene1", "Gene5", "Gene20"),
  BidirectionalSet = data.frame(Gene = c("Gene2", "Gene10", "Gene15"), Direction = c(1, -1, 1))
)

results <- runGSEA(DEGList, gene_sets)
print(results)
#> $Contrast1
#>              pathway      pval      padj    log2err         ES        NES  size
#>               <char>     <num>     <num>      <num>      <num>      <num> <int>
#> 1: UnidirectionalSet 0.6613528 0.6613528 0.02255486 -0.4824193 -0.8464888     3
#> 2:  BidirectionalSet 0.4163603 0.6613528 0.02666464 -0.6122216 -1.0612943     3
#>     leadingEdge stat_used
#>          <list>    <char>
#> 1:        Gene1         B
#> 2: Gene10,Gene2         t
#> 
#> $Contrast2
#>              pathway      pval      padj    log2err        ES       NES  size
#>               <char>     <num>     <num>      <num>     <num>     <num> <int>
#> 1: UnidirectionalSet 0.7828653 0.9450915 0.01824776 0.4192777 0.7422151     3
#> 2:  BidirectionalSet 0.9450915 0.9450915 0.01581049 0.3195876 0.5772290     3
#>            leadingEdge stat_used
#>                 <list>    <char>
#> 1:        Gene1,Gene20         B
#> 2: Gene10,Gene15,Gene2         t
#> 
```
