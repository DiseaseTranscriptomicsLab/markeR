#' Generate a Correlation Heatmap
#'
#' This function computes and visualizes a correlation heatmap for a given set of genes.
#' Optionally, the heatmap can be generated separately for different conditions based on metadata.
#'
#' @param data A numeric matrix or data frame where rows represent genes and columns represent samples.
#' @param metadata Optional. A data frame containing metadata information for the samples.
#'   The first column must contain sample IDs matching the column names in `data`.
#' @param genes A character vector of gene names to include in the heatmap.
#' @param separate.by Optional. A column name from `metadata` to separate heatmaps by different conditions.
#' @param method Character. Correlation method to be used to compute a correlation matrix using the
#'   `corr`function. One of `"pearson"` (default), `"spearman"`, or `"kendall"`.
#' @param plot Logical. If `TRUE` (default), the resulting heatmap(s) will be displayed.
#' @param colorlist Optional. A named list specifying colors for the heatmap gradient.
#'   Default: `list(low = "blue", mid = "white", high = "red")`.
#' @param limits_colorscale Optional. Numeric vector of length 2 specifying limits for the color scale.
#'   Numbers outside this range will be displayed in grey.
#' @param widthTitle Integer. Maximum character width for wrapping plot titles. Default is `16`.
#' @param ncol Optional. Number of columns for arranging multiple plots when `separate.by` is used.
#' @param nrow Optional. Number of rows for arranging multiple plots when `separate.by` is used.
#' @param detailedresults Logical. If `TRUE`, additional details including correlation matrices are returned.
#' @param title Optional. A character string specifying the title of the heatmap.
#'
#' @return A named list with the following elements:
#'   \describe{
#'     \item{`data`}{A data frame containing the transformed correlation values used for plotting.}
#'     \item{`plot`}{A ggplot2 plot object or a combined plot (if `separate.by` is specified).}
#'     \item{`aux`}{A list containing additional data when `detailedresults = TRUE`, including:
#'       \describe{
#'         \item{`method`}{The correlation method used (`"pearson"`, `"spearman"`, or `"kendall"`).}
#'         \item{`corrmatrix`}{A numeric correlation matrix of genes.}
#'         \item{`metadata`}{(Only if `separate.by` is used) A metadata subset for the condition.}
#'         \item{`corrmatrix_ggplot`}{(Only if `separate.by` is used) A melted version of `corrmatrix` in long format for plotting.}
#'         \item{`plot`}{(Only if `separate.by` is used) The individual correlation heatmap for the condition.}
#'       }
#'     }
#'   }
#'
#' @examples
#' \dontrun{
#' data <- matrix(rnorm(100), nrow=10, dimnames = list(paste0("Gene", 1:10), paste0("Sample", 1:10)))
#' metadata <- data.frame(SampleID = paste0("Sample", 1:10), Condition = rep(c("A", "B"), each = 5))
#' CorrelationHeatmap(data, metadata, genes = rownames(data), separate.by = "Condition")
#' }
#'
#' @export
CorrelationHeatmap <- function(data,
                               metadata=NULL,
                               genes,
                               separate.by=NULL,
                               method="pearson",
                               plot=TRUE,
                               colorlist=list(low = "blue", mid = "white", high = "red"),
                               limits_colorscale=NULL,
                               widthTitle=16,
                               ncol=NULL,
                               nrow=NULL,
                               detailedresults=FALSE,
                               title=NULL,
                               cluster_rows=FALSE,
                               cluster_cols=FALSE){

  resultsList <- list()
  resultsList[["data"]]<- list()
  resultsList[["plot"]] <- list()
  resultsList[["aux"]] <- list()

  data <- data[row.names(data) %in% genes,]


  if (!is.null(separate.by) && is.null(metadata)) {

    stop("separate.by is not NULL but metadata is. Please specify metadata.")

  } else if (!is.null(separate.by) && !is.null(metadata)){ # if we have a condition for division

    df_data_merge <- data.frame(NULL)
    resultsList_aux <- list()

    for (condition in unique(metadata[,separate.by])){

      metadata_subset <- subset(metadata, get(separate.by) == condition)
      data_subset <- data[,metadata_subset[,1]] # first column of metadata needs to be sample ID
      data_subset <- log2(data_subset)
      corrmat <- stats::cor(t(data_subset), method=method)
      corrmat_ggplot <- reshape2::melt(corrmat)  # Melt the matrix into long format

      # concatenate everything
      df_data_merge <- as.data.frame(rbind(df_data_merge,
                                           cbind(corrmat_ggplot,condition)))

      plt <- ggplot2::ggplot(corrmat_ggplot, ggplot2::aes(Var1, Var2, fill = value))

      # add tiles
      plt <- plt +
        ggplot2::geom_tile()

      # change colors
      plt <- plt +
        ggplot2::scale_fill_gradient2(low = colorlist[["low"]], high = colorlist[["high"]], mid = colorlist[["mid"]],
                                      midpoint = 0, limits = limits_colorscale, space = "Lab",
                                      name= paste0(method, "'s \ncoefficient"), na.value = "#DBDBDB")





      # change theme
      plt <- plt +
        ggplot2::theme_minimal()

      # Change plot title
      wrapped_title <- wrap_title(condition, width = widthTitle)
      plt <- plt + ggplot2::ggtitle(wrapped_title)

      # change labels
      plt <- plt +
        ggplot2::labs(x = " ", y = " ") +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
                       axis.text.y = ggplot2::element_text(angle = 0, hjust = 1),
                       plot.title = ggplot2::element_text(hjust = 0.5))


      # save to an aux variable that will be used for plotting the combined plot
      resultsList_aux[[condition]] <- plt

      if(detailedresults){

        resultsList[["aux"]][[condition]] <- list()
        resultsList[["aux"]][[condition]][["method"]] <- method
        resultsList[["aux"]][[condition]][["corrmatrix"]] <- corrmat
        resultsList[["aux"]][[condition]][["metadata"]] <- metadata_subset
        resultsList[["aux"]][[condition]][["corrmatrix_ggplot"]]  <- corrmat_ggplot
        resultsList[["aux"]][[condition]][["plot"]] <- plt

      }


    }

    n <- length(resultsList_aux)
    # Determine grid layout
    if (is.null(ncol) && is.null(nrow)) {

      ncol <- ceiling(sqrt(n))
      nrow <- ceiling(n / ncol)

    } else if (is.null(ncol)){

      ncol <- ceiling(n / nrow)

    } else if (is.null(nrow)){

      nrow <- ceiling(n / ncol)

    }

    combined_plot <- ggpubr::ggarrange(plotlist = resultsList_aux, ncol = ncol, nrow = nrow, common.legend = TRUE, align = "h")

    if (!is.null(title)) {

      wrapped_title <- wrap_title(title, width = widthTitle)
      combined_plot <- ggpubr::annotate_figure(combined_plot,
                                               top = grid::textGrob(wrapped_title, gp = grid::gpar(cex = 1.3, fontsize = 10)))
    }

    if(plot) print(combined_plot)

    resultsList[["data"]] <- df_data_merge
    resultsList[["plot"]] <- combined_plot

  } else { #plot everything together if no dividing condition is specified

    data <- log2(data)
    corrmat <- stats::cor(t(data), method=method)
    corrmat_ggplot <- reshape2::melt(corrmat)  # Melt the matrix into long format

    plt <- ggplot2::ggplot(corrmat_ggplot, aes(Var1, Var2, fill = value))

    # add tiles
    plt <- plt +
      ggplot2::geom_tile()

    # change colors
    plt <- plt +
      ggplot2::scale_fill_gradient2(low = colorlist[["low"]], high = colorlist[["high"]], mid = colorlist[["mid"]],
                                    midpoint = 0, limits = limits_colorscale, space = "Lab",
                                    name= paste0(method, "'s \ncoefficient"), na.value = "#DBDBDB")


    # change labels
    plt <- plt +
      ggplot2::labs(x = " ", y = " ") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
                     axis.text.y = ggplot2::element_text(angle = 0, hjust = 1),
                     plot.title = ggplot2::element_text(hjust = 0.5))

    # change theme
    plt <- plt +
      ggplot2::theme_minimal()

    if (!is.null(title)) {

      # Change plot title
      wrapped_title <- wrap_title(title, width = widthTitle)
      plt <- plt + ggplot2::ggtitle(widthTitle)

    }

    # if user wants to print the plot
    if(plot) print(plt)

    resultsList[["data"]] <- corrmat_ggplot
    resultsList[["plot"]] <- plt

    if(detailedresults){

      resultsList[["aux"]][["method"]] <- method
      resultsList[["aux"]][["corrmatrix"]] <- corrmat

    }


  }

  invisible(resultsList)


}
