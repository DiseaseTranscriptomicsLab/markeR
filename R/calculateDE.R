#' Calculate Differential Gene Expression Statistics using limma
#'
#' This function computes differential gene expression statistics for each gene
#' using a linear model via the limma package. Users may supply a custom design
#' matrix directly via the \code{design} argument, or specify a model formula
#' (\code{lmexpression}) (e.g., \code{~0 + X} or \code{~X}) or variables from
#' \code{metadata} to build the design matrix. When contrasts are supplied, they
#' are applied using \code{limma::makeContrasts} and
#' \code{limma::contrasts.fit}. Alternatively, when using \code{lmexpression} or
#' a supplied \code{design}, specific coefficient indices may be provided via
#' \code{coefs} to extract the corresponding gene-level statistics.
#'
#' @param data A numeric matrix of gene expression values with genes as rows and
#'   samples as columns. Row names must correspond to gene identifiers. Data
#'   should *not* be transformed (i.e., not log2 transformed).
#' @param metadata A data frame containing sample metadata used to build the
#'   design matrix (unless a design is provided directly).
#' @param variables A character vector specifying the variable(s) from
#'   \code{metadata} to use in the default linear model. Ignored if
#'   \code{lmexpression} or \code{design} is provided.
#' @param modelmat (Optional) A user-supplied design matrix. If provided, this
#'   design is used directly and \code{lmexpression} and \code{variables} are
#'   ignored. The order of samples in the design matrix should match the order
#'   in data.
#' @param contrasts A character vector specifying contrasts to be applied (e.g.,
#'   \code{c("A-B")}). If multiple contrasts are provided, the function returns
#'   a list of DE results (one per contrast). *Required* if \code{lmexpression}
#'   is NULL, optional otherwise. If not provided, the average expression
#'   profile of each condition will be returned instead of differential gene
#'   expression.
#' @param ignore_NAs Boolean (default: FALSE). Whether to ignore NAs in the
#'   metadata. If TRUE, rows with any NAs will be removed before analysis,
#'   leading to a loss of data to be fitted in the model. Only applicable if
#'   \code{variables} is provided.
#'
#' @return A list of data-frames of differential expression statistics
#'
#' @details The function fits a linear model with \code{limma::lmFit} and
#' applies empirical Bayes moderation with \code{limma::eBayes}. Depending on
#' the input:
#' \itemize{
#'   \item If a design matrix is provided via \code{design}, that design is
#'   used directly.
#'   \item Otherwise, a design matrix is constructed using the \code{variables}
#'   argument (with no intercept).
#'   \item If contrasts are provided, they are applied using
#'   \code{limma::makeContrasts} and \code{limma::contrasts.fit}.
#'   \item If no contrasts are provided, the function returns all possible
#'   coefficients fitted in the linear model.
#' }
#'
#' @examples
#' # Simulate non-negative gene expression data (counts)
#' set.seed(123)
#' expr <- matrix(rpois(1000, lambda = 20), nrow = 100, ncol = 10)
#' rownames(expr) <- paste0("gene", 1:100)
#' colnames(expr) <- paste0("sample", 1:10)
#'
#' # Simulate metadata with a group variable
#' metadata <- data.frame(
#'  sample = colnames(expr),
#'  Group = rep(c("A", "B"), each = 5)
#' )
#'
#' # Differential expression for Group A vs Group B using variables
#' de_var <- calculateDE(
#'   data = expr,
#'   metadata = metadata,
#'   variables = "Group",
#'   contrasts = "A-B"
#' )
#' head(de_var[["A-B"]])
#'
#' # Build equivalent design matrix manually
#' design <- model.matrix(~0 + Group, data = metadata)
#' colnames(design) <- c("A","B")
#'
#' # Differential expression using the design matrix directly
#' de_mat <- calculateDE(
#'   data = expr,
#'   metadata = metadata,
#'   modelmat = design,
#'   contrasts = "A-B"
#' )
#' head(de_mat[["A-B"]])
#'
#' @importFrom limma lmFit eBayes makeContrasts contrasts.fit topTable
#' @export
calculateDE <- function(data, metadata=NULL, variables=NULL, modelmat = NULL,
                        contrasts = NULL, ignore_NAs = FALSE) {

  data <- as.data.frame(data) # Ensure data is a data frame
  remove_prefix <- function(colnames_vector, prefixes) {
    for (prefix in prefixes) {
      new_colnames <- gsub(paste0("^", prefix), "", colnames_vector)

      # Only update values where removal does NOT result in an empty string
      colnames_vector <- ifelse(new_colnames != "", new_colnames, colnames_vector)
    }
    return(colnames_vector)
  }

  # extract variables from equations
  extract_variables <- function(equation) {
    variables <- unique(unlist(strsplit(gsub("[~+*/()]", " ", equation), "\\s+")))
    variables <- setdiff(variables, "")  # Remove empty strings
    return(variables)
  }


  # # Validate inputs
  # if (!is.matrix(data) && !is.data.frame(data)) stop(
  #   "Error: 'data' must be a matrix or a data frame.")
  # if (is.null(rownames(data))) stop(
  #   "Error: 'data' must have row names corresponding to gene identifiers.")
  # if (!is.null(metadata) && !is.data.frame(metadata)) stop(
  #   "Error: 'metadata' must be a data frame.")
  # if (!is.null(metadata) && (ncol(data) != nrow(metadata))) stop(
  #   "Error: Number of samples in 'data' does not match number of rows in 'metadata'.")
  
  # Validate inputs
  if ((!is.matrix(data) && !is.data.frame(data)) || is.null(rownames(data))) {
    stop("Error: 'data' must be a matrix or data frame with row names corresponding to gene identifiers.")
  }
  
  if ((!is.null(metadata) && (!is.data.frame(metadata)) || ncol(data) != nrow(metadata))) {
    stop("Error: 'metadata' must be a data frame with the same number of rows as columns in 'data'.")
  }
  

  # add "." after each variable and remove spaces
  # Important to avoid errors in design matrix
  # if (!is.null(metadata)){
  #   #metadata <- as.data.frame(lapply(metadata, function(x) paste0(".", x)))
  #   #colnames(metadata) <- paste0(colnames(metadata),".")
  #   metadata <- as.data.frame(lapply(metadata, function(x) gsub(" ", "", x)))
  # }
  # if(!is.null(variables)){
  #   variables <- gsub(" ", "", variables)
  #   #variables <- paste0(variables,".")
  # }
  #if(!is.null(modelmat)) colnames(modelmat) <- gsub(" ", "", colnames(modelmat))


  if (!is.null(metadata)){
    # Reorder and subset metadata to match data
    # counts: matrix or data frame with column names as sample IDs
    # metadata: data frame with at least one column containing sample IDs
    
    # 1. Find the metadata column that best matches column names of count matrix
    sample_ids <- colnames(data)
    #best_match_col <- which.max(sapply(metadata, function(col) sum(sample_ids %in% col)))
    best_match_col <- which.max(vapply(metadata, function(col) sum(sample_ids %in% col), numeric(1)))
    
    # print message saying which column was used to match samples 
    message("Using metadata column '", colnames(metadata)[best_match_col], "' to match samples (data column names).")
    
    # 2. Extract matched column
    matched_col <- metadata[[best_match_col]]
    
    # 3. Subset metadata to only those samples present in the count matrix
    metadata_matched <- metadata[matched_col %in% sample_ids, ]
    
    # 4. Reorder metadata to match column order of count matrix
    rownames(metadata_matched) <- metadata_matched[[best_match_col]]
    # drop = FALSE to preserve data frame format
    metadata_matched <- metadata_matched[sample_ids, , drop = FALSE]
    metadata <- metadata_matched
  }


  if (ignore_NAs & !is.null(variables)) {
    # Remove rows with NAs in the specified variables
    na_rows <- apply(metadata[variables], 1, function(x) any(is.na(x)))
    metadata <- metadata[!na_rows, ]
    # Ensure data matches the filtered metadata
    data <- data[, rownames(metadata), drop = FALSE]
  }


  # Construct design matrix
  design_matrix <- tryCatch({

    if (!is.null(modelmat)) {
      # if (!is.matrix(modelmat)) stop("Error: 'modelmat' must be a matrix.")
      # if (nrow(modelmat) != ncol(data)) stop(
      #   "Error: Rows in 'modelmat' must match the number of samples in 'data'.
      #   Check if your metadata has any NAs or consider using ignore_NAs = TRUE.")

      if (!is.matrix(modelmat) || (nrow(modelmat) != ncol(data))) {
        stop("Error: 'modelmat' must be a matrix with rows equal to the number of samples in 'data'.
       Check if your metadata has any NAs or consider using ignore_NAs = TRUE.")
      }

      modelmat
    # } else if (!is.null(lmexpression)) {
    #   lmexpression <- as.formula(lmexpression, env = parent.frame())
    #   design_matrix <- model.matrix(lmexpression, data = metadata)
    #   vars <- extract_variables(lmexpression)
    #   #colnames(design_matrix) <- gsub("^Condition","",colnames(design_matrix))
    #
    #   colnames(design_matrix) <-  remove_prefix(colnames(design_matrix), vars)
    #   colnames(design_matrix) <- gsub(" ", "", colnames(design_matrix))
    #   #colnames(design_matrix) <- sub("^[^.]*\\.", "", colnames(design_matrix))
    #   design_matrix
    } else {
      design_formula <- as.formula(paste("~0+", paste(variables, collapse = " + ")))
      design_matrix <- stats::model.matrix(design_formula, data = metadata)
      #colnames(design_matrix) <- gsub("^Condition","",colnames(design_matrix))
      #colnames(design_matrix) <- sub("^[^.]*\\.", "", colnames(design_matrix)) # remove the variable name
      colnames(design_matrix) <-   remove_prefix(colnames(design_matrix), variables)
      colnames(design_matrix) <- gsub(" ", "", colnames(design_matrix)) # remove spaces
      design_matrix
    }
  }, error = function(e) {
    stop("Error constructing design matrix: ", e$message)
  })


  # Ensure design matrix and data match
  if (nrow(design_matrix) != ncol(data))
    stop("Error: Mismatch between number of samples in 'data' and rows in 'design_matrix'.
         Check if your metadata has any NAs or consider using ignore_NAs = TRUE.")

  # Fit model
  fit <- tryCatch({
    limma::lmFit(log2(data), design_matrix)
  }, error = function(e) {
    stop("Error in lmFit: ", e$message)
  })
  fit <- limma::eBayes(fit)

  resultsList <- list()

  # Process contrasts
  if (!is.null(contrasts)) {

    contrasts <- gsub(" ", "", contrasts) # remove white spaces

    contrast_matrix <- limma::makeContrasts(contrasts = contrasts, levels = design_matrix)
    fit <- limma::contrasts.fit(fit, contrast_matrix)
    fit <- limma::eBayes(fit)

    for (cont in contrasts) {

      if (!(cont %in% colnames(contrast_matrix))) {
        warning("Warning: Contrast '", cont, "' not found in contrast matrix. Skipping.")
        next
      }
      resultsList[[cont]] <-  limma::topTable(fit, cont, n = Inf, sort.by = "logFC")
    }

  } else{

    for (coef in colnames(fit$design)){

      resultsList[[coef]] <-  limma::topTable(fit, coef, n = Inf, sort.by = "logFC")

    }



  }

  return(resultsList)
}
