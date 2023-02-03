

## Perform GWAS with hglm and extra variance components, separating
## the two crossbreds

library(assertthat)
library(dplyr)
library(hglm)
library(purrr)
library(readr)


source("R/hglm_helper_functions.R")

source("R/hglm_gwas_prepare_data_crosses_separate.R")


gwas_pen_bovans_load <- run_gwas(pen_bovans$load_N$data$X6,
                                 model.matrix(~ 1 + weight, pen_bovans_covar$load_N),
                                 cbind(Z_grm_pen_bovans$load_N, Z_pen_bovans_group$load_N),
                                 c(ncol(Z_grm_pen_bovans$load_N), ncol(Z_pen_bovans_group$load_N)),
                                 snps_pruned_pen_bovans$load_N)


gwas_pen_lsl_load <- run_gwas(pen_lsl$load_N$data$X6,
                              model.matrix(~ 1 + weight, pen_lsl_covar$load_N),
                              cbind(Z_grm_pen_lsl$load_N, Z_pen_lsl_group$load_N),
                              c(ncol(Z_grm_pen_lsl$load_N), ncol(Z_pen_lsl_group$load_N)),
                              snps_pruned_pen_lsl$load_N)

saveRDS(gwas_pen_bovans_load,
        file = "gwas/hglm_gwas_pen_bovans_load.Rds")

saveRDS(gwas_pen_lsl_load,
        file = "gwas/hglm_gwas_pen_lsl_load.Rds")



gwas_cage_bovans_load <- run_gwas(cage_bovans$load_N$data$X6,
                          model.matrix(~ 1 + weight, cage_bovans_covar$load_N),
                          Z_grm_cage_bovans$load_N,
                          ncol(Z_grm_cage_bovans$load_N),
                          snps_pruned_cage_bovans$load_N)

gwas_cage_lsl_load <- run_gwas(cage_lsl$load_N$data$X6,
                               model.matrix(~ 1 + weight, cage_lsl_covar$load_N),
                               Z_grm_cage_lsl$load_N,
                               ncol(Z_grm_cage_lsl$load_N),
                               snps_pruned_cage_lsl$load_N)

saveRDS(gwas_cage_bovans_load,
        file = "gwas/hglm_gwas_cage_bovans_load.Rds")

saveRDS(gwas_cage_lsl_load,
        file = "gwas/hglm_gwas_cage_lsl_load.Rds")


gwas_all_bovans_load <- run_gwas(all_bovans$load_N$data$X6,
                                 model.matrix(~ 1 + weight + cage.pen, all_bovans_covar$load_N),
                                 cbind(Z_grm_all_bovans$load_N, Z_all_bovans_group$load_N),
                                 c(ncol(Z_grm_all_bovans$load_N), ncol(Z_all_bovans_group$load_N)),
                                 snps_pruned_all_bovans$load_N)

gwas_all_lsl_load <- run_gwas(all_lsl$load_N$data$X6,
                              model.matrix(~ 1 + weight + cage.pen, all_lsl_covar$load_N),
                              cbind(Z_grm_all_lsl$load_N, Z_all_lsl_group$load_N),
                              c(ncol(Z_grm_all_lsl$load_N), ncol(Z_all_lsl_group$load_N)),
                              snps_pruned_all_lsl$load_N)

saveRDS(gwas_all_bovans_load,
        file = "gwas/hglm_gwas_all_bovans_load.Rds")

saveRDS(gwas_all_lsl_load,
        file = "gwas/hglm_gwas_all_lsl_load.Rds")





gwas_pen_bovans_weight <- run_gwas(pen_bovans$weight$data$X6,
                                 model.matrix(~ 1, pen_bovans_covar$weight),
                                 cbind(Z_grm_pen_bovans$weight, Z_pen_bovans_group$weight),
                                 c(ncol(Z_grm_pen_bovans$weight), ncol(Z_pen_bovans_group$weight)),
                                 snps_pruned_pen_bovans$weight)


gwas_pen_lsl_weight <- run_gwas(pen_lsl$weight$data$X6,
                              model.matrix(~ 1, pen_lsl_covar$weight),
                              cbind(Z_grm_pen_lsl$weight, Z_pen_lsl_group$weight),
                              c(ncol(Z_grm_pen_lsl$weight), ncol(Z_pen_lsl_group$weight)),
                              snps_pruned_pen_lsl$weight)

saveRDS(gwas_pen_bovans_weight,
        file = "gwas/hglm_gwas_pen_bovans_weight.Rds")

saveRDS(gwas_pen_lsl_weight,
        file = "gwas/hglm_gwas_pen_lsl_weight.Rds")



gwas_cage_bovans_weight <- run_gwas(cage_bovans$weight$data$X6,
                                  model.matrix(~ 1, cage_bovans_covar$weight),
                                  Z_grm_cage_bovans$weight,
                                  ncol(Z_grm_cage_bovans$weight),
                                  snps_pruned_cage_bovans$weight)

gwas_cage_lsl_weight <- run_gwas(cage_lsl$weight$data$X6,
                               model.matrix(~ 1, cage_lsl_covar$weight),
                               Z_grm_cage_lsl$weight,
                               ncol(Z_grm_cage_lsl$weight),
                               snps_pruned_cage_lsl$weight)

saveRDS(gwas_cage_bovans_weight,
        file = "gwas/hglm_gwas_cage_bovans_weight.Rds")

saveRDS(gwas_cage_lsl_weight,
        file = "gwas/hglm_gwas_cage_lsl_weight.Rds")


gwas_all_bovans_weight <- run_gwas(all_bovans$weight$data$X6,
                                 model.matrix(~ 1 + cage.pen, all_bovans_covar$weight),
                                 cbind(Z_grm_all_bovans$weight, Z_all_bovans_group$weight),
                                 c(ncol(Z_grm_all_bovans$weight), ncol(Z_all_bovans_group$weight)),
                                 snps_pruned_all_bovans$weight)

gwas_all_lsl_weight <- run_gwas(all_lsl$weight$data$X6,
                              model.matrix(~ 1 + cage.pen, all_lsl_covar$weight),
                              cbind(Z_grm_all_lsl$weight, Z_all_lsl_group$weight),
                              c(ncol(Z_grm_all_lsl$weight), ncol(Z_all_lsl_group$weight)),
                              snps_pruned_all_lsl$weight)

saveRDS(gwas_all_bovans_weight,
        file = "gwas/hglm_gwas_all_bovans_weight.Rds")

saveRDS(gwas_all_lsl_weight,
        file = "gwas/hglm_gwas_all_lsl_weight.Rds")



## pQCT

gwas_pen_bovans_ct <- pmap_dfr(list(data = pen_bovans[3:5],
                                    covar = pen_bovans_covar[3:5],
                                    Z_grm = Z_grm_pen_bovans[3:5],
                                    Z_group = Z_pen_bovans_group[3:5],
                                    snps_pruned = snps_pruned_pen_bovans[3:5]),
                               function(data, covar, Z_grm, Z_group, snps_pruned)
                                       run_gwas(data$data$X6,
                                                model.matrix(~ 1 + weight, covar),
                                                cbind(Z_grm, Z_group),
                                                c(ncol(Z_grm), ncol(Z_group)),
                                                snps_pruned),
                               .id = "trait")


gwas_cage_bovans_ct <- pmap_dfr(list(data = cage_bovans[3:5],
                                     covar = cage_bovans_covar[3:5],
                                     Z_grm = Z_grm_cage_bovans[3:5],
                                     Z_group = Z_cage_bovans_group[3:5],
                                     snps_pruned = snps_pruned_cage_bovans[3:5]),
                                function(data, covar, Z_grm, Z_group, snps_pruned)
                                        run_gwas(data$data$X6,
                                                 model.matrix(~ 1 + weight, covar),
                                                 cbind(Z_grm, Z_group),
                                                 c(ncol(Z_grm), ncol(Z_group)),
                                                 snps_pruned),
                                .id = "trait")


gwas_all_bovans_ct <- pmap_dfr(list(data = all_bovans[3:5],
                                    covar = all_bovans_covar[3:5],
                                    Z_grm = Z_grm_all_bovans[3:5],
                                    Z_group = Z_all_bovans_group[3:5],
                                    snps_pruned = snps_pruned_all_bovans[3:5]),
                               function(data, covar, Z_grm, Z_group, snps_pruned)
                                       run_gwas(data$data$X6,
                                                model.matrix(~ 1 + weight, covar),
                                                cbind(Z_grm, Z_group),
                                                c(ncol(Z_grm), ncol(Z_group)),
                                                snps_pruned),
                               .id = "trait")


saveRDS(gwas_pen_bovans_ct,
        file = "gwas/hglm_gwas_pen_bovans_ct.Rds")

saveRDS(gwas_cage_bovans_ct,
        file = "gwas/hglm_gwas_cage_bovans_ct.Rds")

saveRDS(gwas_all_bovans_ct,
        file = "gwas/hglm_gwas_all_bovans_ct.Rds")



gwas_pen_lsl_ct <- pmap_dfr(list(data = pen_lsl[3:5],
                                 covar = pen_lsl_covar[3:5],
                                 Z_grm = Z_grm_pen_lsl[3:5],
                                 Z_group = Z_pen_lsl_group[3:5],
                                 snps_pruned = snps_pruned_pen_lsl[3:5]),
                            function(data, covar, Z_grm, Z_group, snps_pruned)
                                    run_gwas(data$data$X6,
                                             model.matrix(~ 1 + weight, covar),
                                             cbind(Z_grm, Z_group),
                                             c(ncol(Z_grm), ncol(Z_group)),
                                             snps_pruned),
                            .id = "trait")


gwas_cage_lsl_ct <- pmap_dfr(list(data = cage_lsl[3:5],
                                  covar = cage_lsl_covar[3:5],
                                  Z_grm = Z_grm_cage_lsl[3:5],
                                  Z_group = Z_cage_lsl_group[3:5],
                                  snps_pruned = snps_pruned_cage_lsl[3:5]),
                             function(data, covar, Z_grm, Z_group, snps_pruned)
                                     run_gwas(data$data$X6,
                                              model.matrix(~ 1 + weight, covar),
                                              cbind(Z_grm, Z_group),
                                              c(ncol(Z_grm), ncol(Z_group)),
                                              snps_pruned),
                             .id = "trait")


gwas_all_lsl_ct <- pmap_dfr(list(data = all_lsl[3:5],
                                 covar = all_lsl_covar[3:5],
                                 Z_grm = Z_grm_all_lsl[3:5],
                                 Z_group = Z_all_lsl_group[3:5],
                                 snps_pruned = snps_pruned_all_lsl[3:5]),
                            function(data, covar, Z_grm, Z_group, snps_pruned)
                                    run_gwas(data$data$X6,
                                             model.matrix(~ 1 + weight, covar),
                                             cbind(Z_grm, Z_group),
                                             c(ncol(Z_grm), ncol(Z_group)),
                                             snps_pruned),
                            .id = "trait")


saveRDS(gwas_pen_lsl_ct,
        file = "gwas/hglm_gwas_pen_lsl_ct.Rds")

saveRDS(gwas_cage_lsl_ct,
        file = "gwas/hglm_gwas_cage_lsl_ct.Rds")

saveRDS(gwas_all_lsl_ct,
        file = "gwas/hglm_gwas_all_lsl_ct.Rds")



## Conditional GWAS of major weight locus

peak_marker <- as.character(gwas_all_lsl_weight$marker_id[which.min(gwas_all_lsl_weight$p)])

X_covar <- model.matrix(~ 1 + cage.pen, all_lsl_covar$weight)

X_covar_snp <- cbind(X_covar, as.data.frame(snps_pruned_all_lsl$weight)[, peak_marker])

gwas_conditional <- run_gwas(all_lsl$weight$data$X6,
                             X_covar_snp,
                             cbind(Z_grm_all_lsl$weight, Z_all_lsl_group$weight),
                             c(ncol(Z_grm_all_lsl$weight), ncol(Z_all_lsl_group$weight)),
                             snps_pruned_all_lsl$weight)


saveRDS(gwas_conditional,
        file = "gwas/hglm_gwas_all_lsl_weight_conditional.Rds")
