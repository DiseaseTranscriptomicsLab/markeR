# Metadata for Marthandan et al. (2016) RNA-Seq Study

A data frame containing metadata for samples from the Marthandan et al.
(2016) study (GEO code GSE63577).

## Usage

``` r
data(metadata_example)
```

## Format

A data frame with 45 rows and 6 columns:

- sampleID:

  Unique sample identifier.

- DatasetID:

  Identifier for the dataset (e.g., "Marthandan2016").

- CellType:

  Cell type, e.g. "Fibroblast".

- Condition:

  Experimental condition ("Senescent" or "Proliferative").

- SenescentType:

  Mechanism of senescence (e.g., "Telomere shortening" for senescent
  samples, "none" for proliferative).

- Treatment:

  Treatment or age descriptor (e.g., "PD72 (Replicative senescence)" for
  senescent samples, "young" for proliferative).

## Source

<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63577>

## References

Marthandan S, Priebe S, Baumgart M, Groth M et al. Similarities in Gene
Expression Profiles during In Vitro Aging of Primary Human Embryonic
Lung and Foreskin Fibroblasts. Biomed Res Int 2015;2015:731938. PMID:
26339636

Marthandan S, Baumgart M, Priebe S, Groth M et al. Conserved Senescence
Associated Genes and Pathways in Primary Human Fibroblasts Detected by
RNA-Seq. PLoS One 2016;11(5):e0154531. PMID: 27140416
