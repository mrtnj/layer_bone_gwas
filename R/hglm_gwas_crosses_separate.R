

## Perform GWAS with hglm and extra variance components, separating
## the two crossbreds

library(assertthat)
library(dplyr)
library(hglm)
library(readr)


source("R/hglm_helper_functions.R")

source("R/hglm_gwas_prepare_data_crosses_separate.R")




gwas_pen_bovans_load <- run_gwas(pen_bovans_load$data$X6,
                                 model.matrix(~ 1 + weight, pen_bovans_load_covar),
                                 cbind(Z_grm_pen_bovans_load, Z_pen_bovans_load_group),
                                 c(ncol(Z_grm_pen_bovans_load), ncol(Z_pen_bovans_load_group)),
                                 snps_pruned_pen_bovans_load)


gwas_pen_lsl_load <- run_gwas(pen_lsl_load$data$X6,
                              model.matrix(~ 1 + weight, pen_lsl_load_covar),
                              cbind(Z_grm_pen_lsl_load, Z_pen_lsl_load_group),
                              c(ncol(Z_grm_pen_lsl_load), ncol(Z_pen_lsl_load_group)),
                              snps_pruned_pen_lsl_load)

saveRDS(gwas_pen_bovans_load,
        file = "gwas/hglm_gwas_pen_bovans_load.Rds")

saveRDS(gwas_pen_lsl_load,
        file = "gwas/hglm_gwas_pen_lsl_load.Rds")



gwas_cage_bovans_load <- run_gwas(cage_bovans_load$data$X6,
                          model.matrix(~ 1 + weight, cage_bovans_load_covar),
                          Z_grm_cage_bovans_load,
                          ncol(Z_grm_cage_bovans_load),
                          snps_pruned_cage_bovans_load)

gwas_cage_lsl_load <- run_gwas(cage_lsl_load$data$X6,
                               model.matrix(~ 1 + weight, cage_lsl_load_covar),
                               Z_grm_cage_lsl_load,
                               ncol(Z_grm_cage_lsl_load),
                               snps_pruned_cage_lsl_load)

saveRDS(gwas_cage_bovans_load,
        file = "gwas/hglm_gwas_cage_bovans_load.Rds")

saveRDS(gwas_cage_lsl_load,
        file = "gwas/hglm_gwas_cage_lsl_load.Rds")


gwas_all_bovans_load <- run_gwas(all_bovans_load$data$X6,
                                 model.matrix(~ 1 + weight + cage.pen, all_bovans_load_covar),
                                 cbind(Z_grm_all_bovans_load, Z_all_bovans_load_group),
                                 c(ncol(Z_grm_all_bovans_load), ncol(Z_all_bovans_load_group)),
                                 snps_pruned_all_bovans_load)

gwas_all_lsl_load <- run_gwas(all_lsl_load$data$X6,
                              model.matrix(~ 1 + weight + cage.pen, all_lsl_load_covar),
                              cbind(Z_grm_all_lsl_load, Z_all_lsl_load_group),
                              c(ncol(Z_grm_all_lsl_load), ncol(Z_all_lsl_load_group)),
                              snps_pruned_all_lsl_load)

saveRDS(gwas_all_bovans_load,
        file = "gwas/hglm_gwas_all_bovans_load.Rds")

saveRDS(gwas_all_lsl_load,
        file = "gwas/hglm_gwas_all_lsl_load.Rds")





gwas_pen_bovans_weight <- run_gwas(pen_bovans_weight$data$X6,
                                 model.matrix(~ 1, pen_bovans_weight_covar),
                                 cbind(Z_grm_pen_bovans_weight, Z_pen_bovans_weight_group),
                                 c(ncol(Z_grm_pen_bovans_weight), ncol(Z_pen_bovans_weight_group)),
                                 snps_pruned_pen_bovans_weight)


gwas_pen_lsl_weight <- run_gwas(pen_lsl_weight$data$X6,
                              model.matrix(~ 1, pen_lsl_weight_covar),
                              cbind(Z_grm_pen_lsl_weight, Z_pen_lsl_weight_group),
                              c(ncol(Z_grm_pen_lsl_weight), ncol(Z_pen_lsl_weight_group)),
                              snps_pruned_pen_lsl_weight)

saveRDS(gwas_pen_bovans_weight,
        file = "gwas/hglm_gwas_pen_bovans_weight.Rds")

saveRDS(gwas_pen_lsl_weight,
        file = "gwas/hglm_gwas_pen_lsl_weight.Rds")



gwas_cage_bovans_weight <- run_gwas(cage_bovans_weight$data$X6,
                                  model.matrix(~ 1, cage_bovans_weight_covar),
                                  Z_grm_cage_bovans_weight,
                                  ncol(Z_grm_cage_bovans_weight),
                                  snps_pruned_cage_bovans_weight)

gwas_cage_lsl_weight <- run_gwas(cage_lsl_weight$data$X6,
                               model.matrix(~ 1, cage_lsl_weight_covar),
                               Z_grm_cage_lsl_weight,
                               ncol(Z_grm_cage_lsl_weight),
                               snps_pruned_cage_lsl_weight)

saveRDS(gwas_cage_bovans_weight,
        file = "gwas/hglm_gwas_cage_bovans_weight.Rds")

saveRDS(gwas_cage_lsl_weight,
        file = "gwas/hglm_gwas_cage_lsl_weight.Rds")


gwas_all_bovans_weight <- run_gwas(all_bovans_weight$data$X6,
                                 model.matrix(~ 1 + cage.pen, all_bovans_weight_covar),
                                 cbind(Z_grm_all_bovans_weight, Z_all_bovans_weight_group),
                                 c(ncol(Z_grm_all_bovans_weight), ncol(Z_all_bovans_weight_group)),
                                 snps_pruned_all_bovans_weight)

gwas_all_lsl_weight <- run_gwas(all_lsl_weight$data$X6,
                              model.matrix(~ 1 + cage.pen, all_lsl_weight_covar),
                              cbind(Z_grm_all_lsl_weight, Z_all_lsl_weight_group),
                              c(ncol(Z_grm_all_lsl_weight), ncol(Z_all_lsl_weight_group)),
                              snps_pruned_all_lsl_weight)

saveRDS(gwas_all_bovans_weight,
        file = "gwas/hglm_gwas_all_bovans_weight.Rds")

saveRDS(gwas_all_lsl_weight,
        file = "gwas/hglm_gwas_all_lsl_weight.Rds")

