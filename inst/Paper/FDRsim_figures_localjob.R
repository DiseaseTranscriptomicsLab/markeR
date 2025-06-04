devtools::load_all()


metadata <- readRDS("data/metadata.rds")
corrcounts <- readRDS("data/corrcounts.rds")
signatures_bidirectional <- readRDS("data/SenescenceSigntures_Bidirectional.rds") # Divided by direction

# plt_fdrsim_allsignatures_v2 <- FPR_Simulation(data = corrcounts,
#                                               metadata = metadata,
#                                               original_signatures = signatures_bidirectional,
#                                               gene_list = row.names(corrcounts),
#                                               number_of_sims = 100,
#                                               widthTitle = 30,
#                                               Variable = "Condition",
#                                               titlesize = 12,
#                                               pointSize = 3,
#                                               labsize = 10,
#                                               mode = "simple",
#                                               ColorValues=NULL,
#                                               ncol=NULL,
#                                               nrow=3 )



plt_fdrsim_subset_v2 <- FPR_Simulation(data = corrcounts,
                                       metadata = metadata,
                                       original_signatures = list(CellAge=signatures_bidirectional$CellAge,
                                                                  HernandezSegura=signatures_bidirectional$HernandezSegura,
                                                                  SAUL_SEN_MAYO=signatures_bidirectional$SAUL_SEN_MAYO,
                                                                  SeneQuest=signatures_bidirectional$SeneQuest
                                       ),
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
                                       nrow=1 )


#saveRDS(plt_fdrsim_allsignatures_v2,"data/plt_fdrsim_allsignatures_v2.rds")
saveRDS(plt_fdrsim_subset_v2,"data/plt_fdrsim_subset_v2.rds")
