normalised_counts <- readRDS("data/normalised_counts.rds")
corrcounts <- readRDS("data/corrcounts.rds")
metadata <- readRDS("data/metadata.rds")

library("ggplot2")
library("colorspace")
library("scales")
library("scater")
library("reshape2")

##### Normalised counts

df_metadata <- metadata
mat <- normalised_counts
vars <- c("DatasetID","Condition","CellType")

df_metadata <- df_metadata[, vars, drop = FALSE]

# Ensure mat is genes (rows) x samples (columns)
if (ncol(mat) != nrow(df_metadata)) {
  stop("Number of columns in mat must match number of rows in df_metadata")
}

explained_variance <- sapply(vars, function(var) {
  sapply(1:nrow(mat), function(i) {
    y <- as.numeric(mat[i, ])
    x <- as.factor(df_metadata[[var]])
    model <- lm(y ~ x)
    summary(model)$r.squared
  })
})

# Average variance explained per variable across genes
mean_var <- colMeans(explained_variance, na.rm = TRUE)

df_plot <- data.frame(
  Variable = names(mean_var),
  VarianceExplained = 100 * mean_var
)


saveRDS( df_plot,"data/df_variance_normcounts.rds")



##### Batch corrected counts

df_metadata <- metadata
mat <- corrcounts
vars <- c("DatasetID","Condition","CellType")

df_metadata <- df_metadata[, vars, drop = FALSE]

# Ensure mat is genes (rows) x samples (columns)
if (ncol(mat) != nrow(df_metadata)) {
  stop("Number of columns in mat must match number of rows in df_metadata")
}

explained_variance <- sapply(vars, function(var) {
  sapply(1:nrow(mat), function(i) {
    y <- as.numeric(mat[i, ])
    x <- as.factor(df_metadata[[var]])
    model <- lm(y ~ x)
    summary(model)$r.squared
  })
})

# Average variance explained per variable across genes
mean_var <- colMeans(explained_variance, na.rm = TRUE)

df_plot <- data.frame(
  Variable = names(mean_var),
  VarianceExplained = 100 * mean_var
)


saveRDS( df_plot,"data/df_variance_corrcounts.rds")
