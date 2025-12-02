# Last run: 30 Sept 2025

metadata <- readRDS("metadata.rds")
corrcounts <- readRDS("corrcounts.rds")
signatures_bidirectional <- readRDS("SenescenceSigntures_Bidirectional.rds") # Divided by direction

# Load functions markeR
 
# Path to your folder with functions
fun_dir <- "functions_Sept18_2025"

# Get all .R files in the folder
files <- list.files(paste0("../",fun_dir), pattern = "\\.R$", full.names = TRUE)

# Source them one by one
invisible(lapply(files, source))

 
library("ggplot2")
library("colorspace")
library("scales")
library("scater") 
library("reshape2")
library("ggbreak")
library("ggnewscale")
library("dplyr")
library("tidyr")  
library(patchwork)
library("ggpubr") 
library(purrr)
library(edgeR)
library(RColorBrewer)
library(pROC) 


set.seed("12345678")


FPR_Simulation_adapted <- function(data, metadata, original_signatures, Variable,
                           gene_list = NULL, number_of_sims=100, title=NULL,
                           widthTitle = 30, titlesize = 12,  pointSize = 2,
                           labsize = 10,mode = c( "none","simple","medium","extensive"),
                           ColorValues=NULL, ncol=NULL, nrow=NULL) {
  data <- as.data.frame(data) # Ensure data is a data frame
  if (is.null(gene_list)) gene_list <- row.names(data)
  
  methods <- c("ssGSEA", "logmedian", "ranking")
  
  mode <- match.arg(mode)
  
  type <- identify_variable_type(metadata, Variable)
  
  if(is.null(ColorValues)) ColorValues <- c(Original = "#68B393", Simulated = "#666666")
  
  if (type =="Numeric"){
    
    results <- suppressMessages(CohenF_allConditions(data = data,
                                                     metadata = metadata,
                                                     gene_sets = original_signatures,
                                                     variable = Variable ))
    cohentype <- "f"
    
  } else {
    
    if(mode == "none"){
      
      results <- suppressMessages(CohenF_allConditions(data = data,
                                                       metadata = metadata,
                                                       gene_sets = original_signatures,
                                                       variable = Variable ))
      cohentype <- "f"
      
    } else {
      
      results <- suppressMessages(CohenD_allConditions(data = data,
                                                       metadata = metadata,
                                                       gene_sets = original_signatures,
                                                       variable = Variable, mode = mode))
      cohentype <- "d"
      
    }
  }
  
  
  rows <- list()
  for (signature in names(results)) {
    sig_data <- results[[signature]]
    
    if (cohentype=="d"){
      cohen_mat <- sig_data$CohenD
    } else if (cohentype=="f"){
      cohen_mat <- sig_data$CohenF
    } else {
      stop("Error: results format not valid.")
    }
    
    padj_mat <- sig_data$padj
    
    for (method in rownames(cohen_mat)) {
      for (contrast in colnames(cohen_mat)) {
        rows[[length(rows) + 1]] <- data.frame(
          signature = signature,
          contrast = contrast,
          method = method,
          cohen = abs(cohen_mat[method, contrast]),
          padj = padj_mat[method, contrast],
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  final_df_original_complete <- do.call(rbind, rows)
  
  # Create a list to store plots (one per signature)
  plot_list <- list()
  data_list <- list()
  signature_names <- names(results)
  
  for (sig in signature_names) {
    
    message(paste0("Running simulations for signature "),sig)
    
    # Use the current signature as the base for simulation
    cur_sig <- original_signatures[[sig]]
    
    if (is.vector(cur_sig)) {
      # Treat all genes as "upregulated", will be the same as not giving direction
      cur_sig <- data.frame(Gene = cur_sig, Signal = 1)
    } else {
      colnames(cur_sig) <- c("Gene", "Signal")
    }
    
    
    # Generate simulated signatures based on the current signature
    simulatedsigs <- list()
    #    for (sim in 1:number_of_sims) {
    for (sim in seq_len(number_of_sims)) {
      cur_model_sig <- cur_sig  # copy the current signature
      cur_model_sig$Gene <- sample(gene_list, nrow(cur_sig))  # simulate by sampling genes
      simulatedsigs[[paste0("sim", sim)]] <- cur_model_sig
    }
    
    
    if (cohentype=="d"){
      results2 <- suppressMessages(CohenD_allConditions(
        data = data,
        metadata = metadata,
        gene_sets = simulatedsigs,
        variable = Variable,
        mode = mode
      ))
    } else {
      results2 <- suppressMessages(CohenF_allConditions(
        data = data,
        metadata = metadata,
        gene_sets = simulatedsigs,
        variable = Variable
      ))
    }
    
    
    
    rows <- list()
    for (signature in names(results2)) {
      sig_data <- results2[[signature]]
      
      if (cohentype=="d"){
        cohen_mat <- sig_data$CohenD
      } else if (cohentype=="f"){
        cohen_mat <- sig_data$CohenF
      } else {
        stop("Error: results format not valid.")
      }
      
      padj_mat <- sig_data$padj
      
      for (method in rownames(cohen_mat)) {
        for (contrast in colnames(cohen_mat)) {
          rows[[length(rows) + 1]] <- data.frame(
            signature = signature,
            contrast = contrast,
            method = method,
            cohen = abs(cohen_mat[method, contrast]),
            padj = padj_mat[method, contrast],
            stringsAsFactors = FALSE
          )
        }
      }
    }
    
    final_df_simulated <- do.call(rbind, rows)
    
    
    # Merge Original and Simulated
    
    final_df_original <- final_df_original_complete[final_df_original_complete$signature == sig, ]
    
    final_df_simulated$type <- "Simulated"
    final_df_original$type <- "Original"
    
    final_df <- rbind(final_df_original, final_df_simulated)
    
    # needed to define the quantile dashed lines
    final_df$method <- factor(final_df$method, levels = methods)
    
    
    # Calculate FPR for each Original observation
    final_df$FPR <- NA
    for (i in which(final_df$type == "Original")) {
      row <- final_df[i, ]
      sim_vals <- final_df$cohen[final_df$type == "Simulated" &
                                   final_df$method == row$method &
                                   final_df$contrast == row$contrast]
      fpr <- mean(sim_vals >= row$cohen, na.rm = TRUE)
      final_df$FPR[i] <- fpr
    }
    
    # Define the quantile lines
    methods <- unique(final_df$method)
    contrasts <- unique(final_df$contrast)
    q_data <- data.frame()
    
    for (ct in contrasts) {
      for (mt in methods) {
        subset_df <- final_df[final_df$method == mt & final_df$contrast == ct, ]
        if (nrow(subset_df) == 0) next
        q95 <- stats::quantile(subset_df$cohen, 0.95, na.rm = TRUE)
        xpos <- which(methods == mt)
        q_data <- rbind(q_data, data.frame(
          method = mt, contrast = ct, q_high = q95,
          ypos = xpos, xmin = xpos - 0.3, xmax = xpos + 0.3
        ))
      }
    }
    
    
    # Ensuring the label is always on top
    # Compute max cohen per method + contrast across both Simulated and Original
    all_max <- stats::aggregate(cohen ~ method + contrast, data = final_df, FUN = max)
    
    
    # Extract FPR values from Original rows
    original_df <- final_df[final_df$type == "Original", ]
    
    # Merge FPR values to max data
    all_max$FPR <- NA
    for (i in seq_len(nrow(all_max))) {
      match_idx <- which(original_df$method == all_max$method[i] &
                           original_df$contrast == all_max$contrast[i])
      if (length(match_idx) > 0) {
        all_max$FPR[i] <- original_df$FPR[match_idx[1]]
      }
    }
    
    # Add label text and offset Y position
    all_max$label <- sprintf("FPR=%.2f", all_max$FPR)
    all_max$y <- all_max$cohen + 0.3  # push label above highest point
    
    # Ensure method is factor with same levels as plotting data
    all_max$method <- factor(all_max$method, levels = levels(final_df$method))
    
    # Build the plot for the current signature
    p <- ggplot2::ggplot() +
      geom_jitter(data = final_df[final_df$type == "Simulated",],
                  aes(y = .data$cohen, x = .data$method, color = .data$type),
                  width = 0.3, height = 0,size = pointSize, alpha = 0.5) +
      geom_violin(data = final_df, aes(y = .data$cohen, x = .data$method),
                  fill = "#F0F0F0", color = "black", alpha = 0.5) +
      geom_jitter(data = final_df[final_df$type == "Original",],
                  aes(y = cohen, x = method, color = type),
                  width = 0.3, height = 0, size = pointSize, alpha = 1) +
      geom_text(data = all_max,
                aes(x = .data$method, y = .data$y, label = .data$label),
                size = 3,
                inherit.aes = FALSE) +
      geom_segment(data = q_data,
                   aes(x = .data$xmin, xend = .data$xmax, y = .data$q_high, yend = .data$q_high),
                   linetype = "dashed", color = "red", inherit.aes = FALSE) +
      labs(title = wrap_title(sig, widthTitle),
           y = ifelse(cohentype == "d", "|Cohen's d|", "|Cohen's f|"),
           x = "Method",
           color = "") +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5, size = titlesize),
            axis.text = element_text(size = labsize)) +
      facet_wrap(. ~ contrast, scales = "free", ncol = 1, strip.position = "left") +
      scale_color_manual(values = ColorValues) +
      scale_x_discrete(labels = c(
        logmedian = expression(log[2]*"-median"),
        ranking   = "ranking",
        ssGSEA    = "ssGSEA"
      )) 
    
    plot_list[[sig]] <- p  # ADD THIS LINE to collect each plot
    data_list[[sig]] <- final_df
    
  }  # ADD THIS LINE to close the for-loop over signatures
  
  n <- length(plot_list)
  
  # Determine grid layout.
  if (is.null(ncol) && is.null(nrow)) {
    ncol <- ceiling(sqrt(n))
    nrow <- ceiling(n / ncol)
  } else if (is.null(ncol)) {
    ncol <- ceiling(n / nrow)
  } else if (is.null(nrow)) {
    nrow <- ceiling(n / ncol)
  }
  
  combined_plot <- ggpubr::ggarrange(
    plotlist = plot_list,
    ncol = ncol,
    nrow = nrow,
    common.legend = TRUE,
    legend = "top"
  )
  
  if (!is.null(title)) {
    title <- wrap_title(title, width = widthTitle)
  }
  
  combined_plot <- ggpubr::annotate_figure(
    combined_plot,
    top = grid::textGrob(title, gp = grid::gpar(cex = 1.3, fontsize = titlesize))
  )
  
  print(combined_plot)
  
  invisible(list(plot=combined_plot,
                 data = data_list))
}


plt_fdrsim_subset_test <- FPR_Simulation_adapted(data = corrcounts,
                                         metadata = metadata,
                                         original_signatures = signatures_bidirectional,
                                         gene_list = row.names(corrcounts),
                                         number_of_sims = 100,
                                         widthTitle = 30,
                                         Variable = "Condition",
                                         titlesize = 12,
                                         pointSize = 3,
                                         labsize = 10,
                                         mode = "simple",
                                         ColorValues=NULL,
                                         ncol=NULL,
                                         nrow=3 )


saveRDS(plt_fdrsim_subset_test, "FPRsim_100_data_and_plot_Dec1_2025_allsigs.rds")
