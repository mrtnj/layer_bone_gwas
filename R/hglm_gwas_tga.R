

## Perform GWAS on selected TGA traits

library(assertthat)
library(dplyr)
library(hglm)
library(readr)


source("R/hglm_helper_functions.R")


source("R/hglm_gwas_tga_prepare_data.R")

gwas_pen_load <- run_gwas(pen_load$data$X6,
                          model.matrix(~ 1 + weight + breed, pen_load_covar),
                          cbind(Z_grm_pen_load, Z_pen_load_group),
                          c(ncol(Z_grm_pen_load), ncol(Z_pen_load_group)),
                          snps_pruned_pen_load)

saveRDS(gwas_pen_load,
        file = "gwas/hglm_gwas_pen_load.Rds")


gwas_cage_load <- run_gwas(cage_load$data$X6,
                          model.matrix(~ 1 + weight + breed, cage_load_covar),
                          Z_grm_cage_load,
                          ncol(Z_grm_cage_load),
                          snps_pruned_cage_load)

saveRDS(gwas_cage_load,
        file = "gwas/hglm_gwas_cage_load.Rds")


gwas_all_load <- run_gwas(all_load$data$X6,
                          model.matrix(~ 1 + weight + breed + cage.pen, all_load_covar),
                          cbind(Z_grm_all_load, Z_all_load_group),
                          c(ncol(Z_grm_all_load), ncol(Z_all_load_group)),
                          snps_pruned_all_load)

saveRDS(gwas_all_load,
        file = "gwas/hglm_gwas_all_load.Rds")



gwas_pen_weight <- run_gwas(pen_weight$data$X6,
                            model.matrix(~ 1 + breed, pen_weight_covar),
                            cbind(Z_grm_pen_weight, Z_pen_weight_group),
                            c(ncol(Z_grm_pen_weight), ncol(Z_pen_weight_group)),
                            snps_pruned_pen_weight)

saveRDS(gwas_pen_weight,
        file = "gwas/hglm_gwas_pen_weight.Rds")

gwas_cage_weight <- run_gwas(cage_weight$data$X6,
                             model.matrix(~ 1 + breed, cage_weight_covar),
                             Z_grm_cage_weight,
                             ncol(Z_grm_cage_weight),
                             snps_pruned_cage_weight)

saveRDS(gwas_cage_weight,
        file = "gwas/hglm_gwas_cage_weight.Rds")

gwas_all_weight <- run_gwas(all_weight$data$X6,
                            model.matrix(~ 1 + breed + cage.pen, all_weight_covar),
                            cbind(Z_grm_all_weight, Z_all_weight_group),
                            c(ncol(Z_grm_all_weight), ncol(Z_all_weight_group)),
                            snps_pruned_all_weight)

saveRDS(gwas_all_weight,
        file = "gwas/hglm_gwas_all_weight.Rds")



peak_marker <- as.character(gwas_all_weight$marker_id[which.min(gwas_all_weight$p)])

X_covar <- model.matrix(~ 1 + breed + cage.pen, all_weight_covar)

X_covar_snp <- cbind(X_covar, as.data.frame(snps_pruned_all_weight)[, peak_marker])

gwas_all_weight_conditional <- run_gwas(all_weight$data$X6,
                                        X_covar_snp,
                                        cbind(Z_grm_all_weight, Z_all_weight_group),
                                        c(ncol(Z_grm_all_weight), ncol(Z_all_weight_group)),
                                        snps_pruned_all_weight)

saveRDS(gwas_all_weight_conditional,
        file = "gwas/hglm_gwas_all_weight_conditional.Rds")
