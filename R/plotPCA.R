#' Principal Component Analysis (PCA) Plot
#'
#' This function performs PCA on a given dataset and visualizes the results using ggplot2.
#' It allows users to specify genes of interest, customize scaling and centering, and color points
#' based on a metadata variable.
#'
#' @param data A numeric matrix or data frame where rows represent genes and columns represent samples.
#' @param metadata A data frame containing sample metadata. The first column should contain sample names. Default is NULL.
#' @param genes A character vector specifying genes to be included in the PCA. Default is NULL (uses all genes).
#' @param scale Logical; if TRUE, variables are scaled before PCA. Default is FALSE.
#' @param center Logical; if TRUE, variables are centered before PCA. Default is TRUE.
#' @param PCs A list specifying which principal components (PCs) to plot. Default is \code{list(c(1,2))}.
#' @param ColorVariable A character string specifying the metadata column used for coloring points. Default is NULL.
#' @param ColorValues A vector specifying custom colors for groups in \code{ColorVariable}. Default is NULL.
#' @param pointSize Numeric; sets the size of points in the plot. Default is 5.
#' @param legend_nrow Integer; number of rows in the legend. Default is 2.
#' @param legend_position Character; position of the legend ("bottom", "top", "right", "left"). Default is "bottom".
#' @param ncol Integer; number of columns in the arranged PCA plots. Default is determined automatically.
#' @param nrow Integer; number of rows in the arranged PCA plots. Default is determined automatically.
#'
#' @return An invisible list containing:
#' \describe{
#'   \item{\code{plt}}{A ggplot2 or ggarrange object displaying the PCA plot.}
#'   \item{\code{data}}{A data frame containing PCA-transformed data and sample metadata (if not NULL).}
#' }
#'
#' @details
#' The function performs PCA using \code{prcomp()} and visualizes the results using \code{ggplot2}.
#' If a metadata data frame is provided, it ensures the sample order matches between data and metadata.
#'
#' @examples
#' \dontrun{
#' # Example dataset
#' set.seed(123)
#' data <- matrix(rnorm(1000), nrow=50, ncol=20)
#' colnames(data) <- paste0("Sample", 1:20)
#' rownames(data) <- paste0("Gene", 1:50)
#'
#' metadata <- data.frame(Sample = colnames(data),
#'                        Group = rep(c("A", "B"), each = 10))
#'
#' # Basic PCA plot
#' plotPCA(data, metadata, ColorVariable = "Group")
#' }
#'
#' @importFrom edgeR DGEList
#' @importFrom stats prcomp
#' @import ggplot2
#' @importFrom ggpubr ggarrange
#' @export
#'
plotPCA <- function(data, metadata=NULL, genes=NULL, scale=FALSE, center=TRUE, PCs=list(c(1,2)), ColorVariable=NULL,ColorValues=NULL,pointSize=5,legend_nrow=2, legend_position=c("bottom","top","right","left"),ncol=NULL, nrow=NULL){

  legend_position <- match.arg(legend_position)

  if (is.null(metadata) & !is.null(ColorVariable)) stop("ColorVariable only available when metadata is specified.")
  if (is.null(genes)){
    genes <-  row.names(data)
  }

  data <- data[row.names(data) %in% genes, , drop=F]

  if (!nrow(data)>1) stop(paste0("Error: Number of genes should be >1; In your data you have only found the gene ",genes))

  # Ensure metadata matches sample order if provided
  if (!is.null(metadata)) {
    colnames(metadata)[1] <- "Sample"
    rownames(metadata) <- metadata$Sample
    metadata <- metadata[colnames(data), , drop = FALSE]
    y <- edgeR::DGEList( log2(data+1), samples= metadata)
  } else {
    y <- edgeR::DGEList( log2(data+1))
  }



  nPCs <- max(unlist(PCs)) # get the maximum number of PC based on the user's choice

  PCAdata <- stats::prcomp(t(y$counts), scale=scale, center=center)
  PCAcounts <- PCAdata$x
  PCAcounts <- as.data.frame(PCAcounts)

  if (nPCs > ncol(PCAcounts)) stop("Error: Number of genes too low for number of chosen PCs. Please reduce number of PCs.")

  PCAcounts <-  cbind(PCAcounts[,1:nPCs],y$samples)

  pltList <- list()

  for (pc in PCs){
    pc <- unlist(pc)
    ev = PCAdata$sdev^2
    pc_x <- round(100*ev[pc[1]]/sum(ev),2)
    pc_y <- round(100*ev[pc[2]]/sum(ev),2)

    plt <- ggplot2::ggplot(PCAcounts, ggplot2::aes_string(y = paste0("PC",pc[1]), x =  paste0("PC",pc[2])))


    # Add jittered points, optionally colored by ColorVariable.Default: Brewer Pallette "Paired"
    if (!is.null(ColorVariable)) {
      plt <- plt + ggplot2::geom_point(ggplot2::aes_string(fill = ColorVariable), size = pointSize, alpha = 0.5, shape=21, color="black")
    } else {
      plt <- plt + ggplot2::geom_point(size = pointSize, alpha = 0.5, shape=21, color="black", fill="#D8D8D8")
    }

    # If ColorValues is provided, use a manual color scale; otherwise, if ColorVariable is provided,
    # use a default brewer palette.
    if (!is.null(ColorValues)) {
      plt <- plt + ggplot2::scale_fill_manual(values = ColorValues)
    } else if (!is.null(ColorVariable)) {
      plt <- plt + ggplot2::scale_fill_brewer(palette = "Paired")
    }

    # Add axis labels (including variance)

    xlab <- paste0("PC",pc[1],": ",pc_x,"% variance")
    ylab <- paste0("PC",pc[2],": ",pc_y,"% variance")
    titleplot <- paste0("PC",pc[1], " vs PC",pc[2])

    plt <- plt + ggplot2::labs(fill = "", x = xlab, y = ylab, title=titleplot)

    # Change theme
    plt <- plt +
      ggplot2::theme_bw()+
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust=1),
                     plot.title = ggplot2::element_text(hjust = 0.5),
                     legend.position = "bottom")

    # Adjust legend rows if legend_nrow is specified
    if (!is.null(legend_nrow)) {
      plt <- plt + ggplot2::guides(fill = ggplot2::guide_legend(nrow = legend_nrow, position = legend_position))
    }

    # Add reference lines
    plt <- plt +
      ggplot2::geom_vline(xintercept=0, linetype="dotted") +
      ggplot2::geom_hline(yintercept=0, linetype="dotted")


    pltList <- c(pltList, list(plt))



  }

  n <- length(pltList)

  if(n==1){
    plt <- pltList[[1]]
  } else {

    # Determine grid layout
    if (is.null(ncol) && is.null(nrow)) {

      ncol <- ceiling(sqrt(n))
      nrow <- ceiling(n / ncol)

    } else if (is.null(ncol)){

      ncol <- ceiling(n / nrow)

    } else if (is.null(nrow)){

      nrow <- ceiling(n / ncol)

    }

    if (!is.null(ColorVariable)) {
      plt <- ggpubr::ggarrange(plotlist = pltList, ncol = ncol, nrow = nrow, common.legend = TRUE, align = "h")
    } else {
      plt <- ggpubr::ggarrange(plotlist = pltList, ncol = ncol, nrow = nrow, align = "h")
    }



  }

  print(plt)

  invisible(list(plt = plt,
                 data = PCAcounts))

}
