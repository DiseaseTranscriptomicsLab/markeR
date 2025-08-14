
<!-- README.md is generated from README.Rmd. Please edit that file -->

# markeR

This repository contains all scripts, data references, and resources to
reproduce the analyses and figures presented in the manuscript exploring
senescence-associated gene sets using the **`markeR`** package.

------------------------------------------------------------------------

## Repository Structure

The repository is organized into two main folders: **`preprocessing`**
and **`Figures`**.

- **`preprocessing`**: Contains scripts to download, align,
  quality-check, filter, normalize, and batch-correct RNA-seq datasets,
  producing the processed data used for analysis.
  - `Download_Alignment_QC`: Includes `Metadata.csv`,
    `pipeline_download_processing.sh`, and `sampleIDs.txt` for
    retrieving raw datasets from GEO and ArrayExpress and performing
    alignment and quality control.  
  - `Filtering_Normalisation_Batch`: Includes
    `Preprocessing_SenescenceDatasets.Rmd`, which performs
    low-expression filtering, normalization, and batch effect
    correction.
- **`Figures`**: Contains scripts and outputs for generating manuscript
  figures.
  - `Figures_Paper_Main.Rmd` and `Figures_Paper_Supplementary.Rmd`
    contain the figure-generation code.  
  - `Figs/` stores the resulting `.png` images.

------------------------------------------------------------------------

## Table of Contents

1.  [Data Availability](#data-availability)  
2.  [Alignment](#alignment)  
3.  [Pre-processing](#pre-processing)  
4.  [Generating Figures](#generating-figures)  
5.  [Installation](#installation)  
6.  [References](#references)

------------------------------------------------------------------------

## Data Availability

No data has been generated for this manuscript. The senescence
compendium was compiled using publicly available RNA-seq datasets
accessed under the following accession codes:
`E-MTAB-5403, E-MTAB-9714, GSE175533, GSE60340, GSE63577, GSE64553, GSE214410, GSE222400, GSE224070, GSE230181, GSE235768, GSE247831, GSE74324, GSE75643, GSE94928, GSE130727, GSE250224, GSE196724, GSE206677, GSE217718, GSE110268, 200213323, E-MTAB-10969`.
This data includes six cell types (fibroblast \[n=384\], endothelial
\[n=90\], keratinocyte \[n=36\], melanocyte \[n=12\], mesenchymal
\[n=12\], and neuronal \[n=11\]) and three cellular states: senescent (n
= 240), proliferative (n = 272), and quiescent (n = 33).

Data from the Genotype-Tissue Expression (GTEx) project v8, already
aligned, filtered, normalized, and batch-corrected, was obtained from
the repository associated with Schneider et al., 2024 ([voyAGEr
GitHub](https://github.com/DiseaseTranscriptomicsLab/voyAGEr)). Donor
metadata, including exact age, was retrieved via dbGaP accession number
`phs000424.v9.p2`.

Senescence gene sets that reflect both community use and diversity of
biological contexts were selected, based publicly available resources
available as of 13 March 2024. We included five well-cited, open-source
gene sets derived from different cell types and senescence-inducing
stimuli: `SAUL_SEN_MAYO` (Saul et al., 2022), `CSGene` (Zhao et al.,
2016), `CellAge` (Chatsirisupachai et al., 2019), `SeneQuest` (Gorgoulis
et al., 2019) and `HernandezSegura` (Hernandez-Segura et al., 2017).
Additionally, we incorporated four gene sets from MSigDB (Subramanian et
al., 2005): `GOBP_CELLULAR_SENESCENCE`,
`GOBP_NEGATIVE_REGULATION_OF_CELLULAR_SENESCENCE`,
`GOBP_POSITIVE_REGULATION_OF_CELLULAR_SENESCENCE` and
`REACTOME_CELLULAR_SENESCENCE`.

------------------------------------------------------------------------

## Alignment

Raw FASTQ files from GEO were downloaded using `fasterq-dump` and
manually from ArrayExpress. Quality control was performed using `fastqc`
and summarized with `MultiQC`. All reads were pseudo-aligned to the
RefSeq human transcriptome (NCBI release 109) using `kallisto` v0.44.0.

------------------------------------------------------------------------

## Pre-processing

Lowly expressed genes were removed by retaining only genes with an
average of \>70 reads in at least one condition (proliferative,
quiescent, or senescent). Normalization factors were calculated using
`calcNormFactors` from **edgeR**, and counts were log2-transformed with
**voom** from **limma**. Batch effects were corrected using a modified
`removeBatchEffect` function from **limma**, regressing out dataset ID
while preserving biological variation (e.g., cell type).

The **`preprocessing`** folder contains the full scripts used to
reproduce these steps:  
- `pipeline_download_processing.sh` in the **`Download_Alignment_QC`**
folder automates FASTQ download and alignment.  
- `Preprocessing_SenescenceDatasets.Rmd` in the
**`Filtering_Normalisation_Batch`** performs filtering, normalization,
and batch correction.

------------------------------------------------------------------------

## Generating Figures

Figures for the manuscript can be generated using the R Markdown scripts
in the **`Figures`** folder:  
- `Figures_Paper_Main.Rmd` generates all main manuscript figures.  
- `Figures_Paper_Supplementary.Rmd` generates supplementary figures.

The resulting images are stored in the `Figs` subfolder. These scripts
rely on the processed datasets generated in the **preprocessing**
workflow, not included in this repository.

------------------------------------------------------------------------

## Contact

📩 For any questions or concerns, feel free to reach out:

**Rita Martins-Silva**  
Email: <rita.silva@medicina.ulisboa.pt>
