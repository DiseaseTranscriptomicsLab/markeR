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
#' A numeric matrix containing gene expression data from the Marthandan et al. (2016) study (GEO code GSE63577).
#' The data represent normalized (non-transformed) counts for genes across samples.
#'
#' The samples were quality-controlled using FastQC and aligned using kallisto (v0.46.1) against the
#' RefSeq v109 transcriptome. The original dataset had 45 samples, but intermediate time points from
#' HFF and MRC5 were removed.
#'
#' @format A numeric matrix with rows corresponding to genes and columns corresponding to samples.
#'   Row names indicate gene symbols and column names are sample identifiers.
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63577}
#'
#' @references  Marthandan S, Priebe S, Baumgart M, Groth M et al. Similarities in Gene Expression Profiles during
#' In Vitro Aging of Primary Human Embryonic Lung and Foreskin Fibroblasts. Biomed Res Int 2015;2015:731938. PMID: 26339636
#' @references  Marthandan S, Baumgart M, Priebe S, Groth M et al. Conserved Senescence Associated Genes and Pathways
#' in Primary Human Fibroblasts Detected by RNA-Seq. PLoS One 2016;11(5):e0154531. PMID: 27140416
#'
#' @keywords datasets
"counts_example"
