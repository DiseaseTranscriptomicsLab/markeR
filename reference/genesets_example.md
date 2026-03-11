# Example Gene Sets for Cellular Senescence

Example Gene Sets for Cellular Senescence

## Usage

``` r
data(genesets_example)
```

## Format

A named list of length 3:

- Literature_Senescence:

  Character vector of gene symbols. A small, curated gene set of
  commonly reported senescence markers, with directionality (+1 or -1).

- REACTOME_Senescence:

  Character vector of gene symbols. The REACTOME_CELLULAR_SENESCENCE
  from MSigDB database No directionality.

- HernandezSegura:

  A data frame with columns `gene` and `direction`. A gene set from
  Hernandez-Segura et al. (2017), with directionality (+1 or -1).

## References

Hernandez-Segura A, de Jong TV, Melov S, Guryev V, Campisi J, Demaria M.
Unmasking Transcriptional Heterogeneity in Senescent Cells. *Curr Biol.*
2017 Sep 11;27(17):2652-2660.e4. doi: 10.1016/j.cub.2017.07.033. Epub
2017 Aug 30. PMID: 28844647; PMCID: PMC5788810.
