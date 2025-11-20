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

plt_fdrsim_subset_test <- FPR_Simulation(data = corrcounts,
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


saveRDS(plt_fdrsim_subset_test, "FPRsim_100_data_and_plot_Sept30_2025_allsigs.rds")
