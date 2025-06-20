#' Metadata for Marthandan et al. (2016) RNA-Seq Study
#'
#' A data frame containing metadata for samples from the Marthandan et al. (2016) study (GEO code GSE63577).
#'
#' @format A data frame with 45 rows and 6 columns:
#' \describe{
#'   \item{sampleID}{Unique sample identifier.}
#'   \item{DatasetID}{Identifier for the dataset (e.g., "Marthandan2016").}
#'   \item{CellType}{Cell type, e.g. "Fibroblast".}
#'   \item{Condition}{Experimental condition ("Senescent" or "Proliferative").}
#'   \item{SenescentType}{Mechanism of senescence (e.g., "Telomere shortening" for senescent samples, "none" for proliferative).}
#'   \item{Treatment}{Treatment or age descriptor (e.g., "PD72 (Replicative senescence)" for senescent samples, "young" for proliferative).}
#' }
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63577}
#'
#' @references  Marthandan S, Priebe S, Baumgart M, Groth M et al. Similarities in Gene Expression Profiles during In Vitro Aging of
#' Primary Human Embryonic Lung and Foreskin Fibroblasts. Biomed Res Int 2015;2015:731938. PMID: 26339636
#' @references  Marthandan S, Baumgart M, Priebe S, Groth M et al. Conserved Senescence Associated Genes and Pathways in Primary
#' Human Fibroblasts Detected by RNA-Seq. PLoS One 2016;11(5):e0154531. PMID: 27140416
#'
#' @keywords datasets
"metadata_example"

#' Gene Expression Counts for Marthandan et al. (2016) RNA-Seq Data
#'
#' A numeric matrix containing filtered and normalized gene expression data from the Marthandan et al. (2016) study (GEO accession GSE63577).
#'
#' Raw FASTQ files were downloaded using `fasterq-dump` (v2.11.0) and processed in a reproducible conda environment (Python v3.11.5).
#' Quality control was conducted using FastQC (v0.12.1) and summarised with MultiQC (v1.14). Pseudo-alignment to the RefSeq transcriptome (NCBI release 109)
#' was performed using kallisto (v0.44.0). Genes with low expression (mean count < 70 in all conditions) were filtered out. Count normalization factors were
#' calculated with `edgeR::calcNormFactors`, and log2-transformed values were obtained via `limma::voom`.
#'
#' Intermediate time points for HFF and MRC5 cell lines were excluded, resulting in a final dataset with 45 high-quality samples across proliferative, quiescent,
#' and senescent conditions.
#'
#' @format A numeric matrix with rows as genes (gene symbols) and columns as samples (sample IDs).
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63577}
#'
#' @references Marthandan S, Priebe S, Baumgart M, Groth M et al. Similarities in Gene Expression Profiles during In Vitro Aging of Primary Human Embryonic Lung and
#' Foreskin Fibroblasts. *Biomed Res Int* 2015;2015:731938. PMID: 26339636
#' @references Marthandan S, Baumgart M, Priebe S, Groth M et al. Conserved Senescence Associated Genes and Pathways in Primary Human Fibroblasts Detected by RNA-Seq.
#' *PLoS One* 2016;11(5):e0154531. PMID: 27140416
#'
#' @keywords datasets
"counts_example"


#' Example Gene Sets for Cellular Senescence
#'
#' @format A named list of length 3:
#' \describe{
#'   \item{Literature_Senescence}{Character vector of gene symbols. A small, curated gene set of commonly reported senescence markers, with directionality (+1 or -1).}
#'   \item{REACTOME_Senescence}{Character vector of gene symbols. A gene set from the REACTOME_CELLULAR_SENESCENCE MSigDB pathway. No directionality.}
#'   \item{HernandezSegura}{A data frame with columns `gene` and `direction`. A gene set from Hernandez-Segura et al. (2017), with directionality (+1 or -1).}
#' }
#'
#' @usage data(genesets_example)
#' @keywords datasets
#' @references Hernandez-Segura A, de Jong TV, Melov S, Guryev V, Campisi J, Demaria M. Unmasking Transcriptional Heterogeneity in Senescent Cells.
#' *Curr Biol.* 2017 Sep 11;27(17):2652-2660.e4. doi: 10.1016/j.cub.2017.07.033. Epub 2017 Aug 30. PMID: 28844647; PMCID: PMC5788810.
"genesets_example"

