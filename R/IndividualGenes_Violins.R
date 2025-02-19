#' Generate Violin Plots for Individual Genes
#'
#' This function creates violin plots of gene expression data with jittered points and optional faceting.
#' It allows visualization of individual gene expression distributions across sample groups.
#'
#' @param data A data frame containing gene expression values with row names as gene names and column names as sample IDs. **(Required)**
#' @param metadata An optional data frame containing sample metadata. The first column must match the sample IDs from `data`. **(Optional)**
#' @param genes A character vector of gene names to be plotted.  **(Required)**
#' @param GroupingVariable A character string specifying the column in `metadata` used for grouping samples on the x-axis. **(Required)**
#' @param plot A logical value indicating whether to print the plot. If `FALSE`, only the output list is returned. Default is `TRUE`. **(Optional)**
#' @param ncol An optional numeric value specifying the number of columns in the facet grid. If not provided, it is computed automatically. Only applicable if `divide` is `NULL`. **(Optional)**
#' @param nrow An optional numeric value specifying the number of rows in the facet grid. If not provided, it is computed automatically. Only applicable if `divide` is `NULL`. **(Optional)**
#' @param divide An optional character string specifying a column in `metadata` to be used for facetting, besides faceting by genes. **(Optional)**
#' @param invert_divide A logical value indicating whether to invert the facet layout, when `divide` is being used. Default is `FALSE`, corresponding to genes in the rows. **(Optional)**
#' @param ColorValues An optional named vector mapping unique values of `ColorVariable` to specific colors. If `NULL`, a default Brewer palette ("Paired") is used. **(Optional)**
#' @param pointSize A numeric value specifying the size of the points in the plot. Default is `2`. **(Optional)**
#' @param ColorVariable A character string specifying a metadata column used for coloring points. Default is `NULL`. **(Optional)**
#' @param title A character string specifying the title of the plot. Default is `NULL`. **(Optional)**
#' @param widthTitle A numeric value specifying the maximum width of the title before inserting line breaks. **(Optional)**
#' @param y_limits A numeric vector of length 2 specifying the limits of the y-axis. If `NULL` (default), the y-axis is adjusted automatically. **(Optional)**
#' @param legend_nrow A numeric value specifying the number of rows in the legend. Default is `NULL`. **(Optional)**
#' @param xlab A character string specifying the x-axis label. If `NULL`, it defaults to `GroupingVariable`. **(Optional)**
#' @param colorlab A character string specifying the legend title for colors. Default is an empty string. **(Optional)**
#'
#' @return A list containing:
#'   \item{plot}{A ggplot2 object representing the facetted violin plots.}
#'   \item{data}{A data frame used for plotting, including transformed expression values (log2) and metadata.}
#'
#' @details
#' The function processes the gene expression data, filters for the specified genes, and transforms expression values using `log2()`.
#' A violin plot with jittered points is generated using `ggplot2`. A median summary is added as a crossbar. If `divide` is provided,
#' facets are created using `ggh4x::facet_grid2()`.
#' Color customization is available via `ColorVariable` and `ColorValues`.
#'
#' @export
#' @examples
#' # Example dataset
#' data <- data.frame(
#'   A = c(10, 20, 30),
#'   B = c(5, 15, 25),
#'   C = c(2, 12, 22)
#' )
#' rownames(data) <- c("Gene1", "Gene2", "Gene3")
#'
#' metadata <- data.frame(
#'   sample = c("A", "B", "C"),
#'   Group = c("Control", "Treatment", "Control")
#' )
#'
#' genes <- c("Gene1", "Gene2")
#'
#' IndividualGenes_Violins(data, metadata, genes, "Group")
#'
IndividualGenes_Violins <- function(data, metadata=NULL, genes,GroupingVariable, plot=TRUE, ncol=NULL, nrow=NULL, divide=NULL, invert_divide=FALSE, ColorValues=NULL, pointSize=2, ColorVariable=NULL,title=NULL, widthTitle=16,y_limits = NULL, legend_nrow = NULL,xlab=NULL, colorlab=NULL){

  set.seed("1234")

  # Wrap the signature name using the helper function
  wrapped_title <- wrap_title(title, width = widthTitle)

  datasub <- data[row.names(data) %in% genes,]
  n <- nrow(data[row.names(data) %in% genes,])

  datasub$gene <- row.names(datasub)
  datasub <- reshape2::melt((datasub))
  datasub$value <- log2(datasub$value)
  colnames(datasub) <- c("gene","sample","expression")


  if (!is.null(metadata)) colnames(metadata)[1] <- "sample"
  if(!is.null(metadata)) datasub <- merge(datasub,metadata, by="sample")

  # Determine grid layout
  if (is.null(ncol) && is.null(nrow)) {

    ncol <- ceiling(sqrt(n))
    nrow <- ceiling(n / ncol)

  } else if (is.null(ncol)){

    ncol <- ceiling(n / nrow)

  } else if (is.null(nrow)){

    nrow <- ceiling(n / ncol)

  }

  plt <- ggplot2::ggplot(datasub, ggplot2::aes_string(y = "expression",x=GroupingVariable))

  # Add jittered points, optionally colored by ColorVariable.Default: Brewer Pallette "Paired"
  if (!is.null(ColorVariable)) {
    plt <- plt + ggplot2::geom_jitter(ggplot2::aes_string(color = ColorVariable), size = pointSize, alpha = 0.5)
  } else {
    plt <- plt + ggplot2::geom_jitter(size = pointSize, alpha = 0.5) + ggplot2::scale_color_brewer(palette = "Paired")
  }


  # If ColorValues is provided, use a manual color scale; otherwise, if ColorVariable is provided,
  # use a default brewer palette.
  if (!is.null(ColorValues)) {
    plt <- plt + ggplot2::scale_color_manual(values = ColorValues)
  } else if (!is.null(ColorVariable)) {
    plt <- plt + ggplot2::scale_color_brewer(palette = "Paired")
  }

  # Adjust legend rows if legend_nrow is specified
  if (!is.null(legend_nrow)) {
    plt <- plt + ggplot2::guides(color = ggplot2::guide_legend(nrow = legend_nrow))
  }

  # Add violin
  plt <- plt + ggplot2::geom_violin(alpha=0.4)

  # Add median summary crossbar.
  plt <- plt + ggplot2::stat_summary(fun = median, fun.min = median, fun.max = median,
                                 geom = "crossbar", width = 0.25,
                                 position = ggplot2::position_dodge(width = 0.13))

  # add second layer of facet, if any
  if (!is.null(divide)){
    if (!invert_divide){
      plt <- plt +
        ggh4x::facet_grid2( gene ~ .data[[divide]],
                            scales = "free",
                            independent= "y")
    } else {
      plt <- plt +
        ggh4x::facet_grid2( .data[[divide]] ~ gene,
                            scales = "free",
                            independent= "y") # labeller = label_wrap_gen(width = 2, multi_line = TRUE)
    }

  } else{

    plt <- plt +
      ggplot2::facet_wrap( . ~ gene,
                           scales = "free",
                           ncol=ncol,
                           nrow=nrow)

  }

  # If y_limits is specified, crop the plot without adjusting the data (violins).

  if (!is.null(y_limits)) {
    plt <- plt + ggplot2::coord_cartesian(ylim = y_limits)
  }

  # Change labels
  if(is.null(xlab)){
    xlab <- GroupingVariable
  } else if (is.null(colorlab)){
    colorlab <- ""
  }

  plt <- plt + ggplot2::labs(title = title, color = colorlab, x = xlab, y = "Expression (log2CPM)")

  # change theme
  plt <- plt + ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust=1),
          plot.title = ggplot2::element_text(hjust = 0.5),
          legend.position = "bottom")

  if(plot){
    print(plt)
  }

  invisible(list(plot=plt,
                 data=datasub))
}
