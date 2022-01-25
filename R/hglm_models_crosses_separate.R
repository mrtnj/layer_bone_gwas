

## Genomic models with hglm for separated crosses

library(assertthat)
library(dplyr)
library(hglm)
library(readr)

source("R/hglm_helper_functions.R")


source("R/hglm_gwas_prepare_data_crosses_separate.R")



pen_bovans_load_model <- hglm(y = pen_bovans_load$data$X6,
                              X = model.matrix(~ 1 + weight, pen_bovans_load_covar),
                              Z = cbind(Z_grm_pen_bovans_load, Z_pen_bovans_load_group),
                              RandC = c(ncol(Z_grm_pen_bovans_load), ncol(Z_pen_bovans_load_group)),
                              calc.like = TRUE)

pen_bovans_load_null <- hglm(y = pen_bovans_load$data$X6,
                             X = model.matrix(~ 1 + weight, pen_bovans_load_covar),
                             Z = Z_pen_bovans_load_group,
                             calc.like = TRUE)

lrt(pen_bovans_load_model, pen_bovans_load_null)

cage_bovans_load_model <- hglm(y = cage_bovans_load$data$X6,
                               X = model.matrix(~ 1 + weight, cage_bovans_load_covar),
                               Z = Z_grm_cage_bovans_load,
                               calc.like = TRUE)


lrt(cage_bovans_load_model)

pen_bovans_load_model$SummVC2
get_h2(pen_bovans_load_model)
get_group_ratio(pen_bovans_load_model)

cage_bovans_load_model$SummVC2
get_h2(cage_bovans_load_model)



pen_lsl_load_model <- hglm(y = pen_lsl_load$data$X6,
                              X = model.matrix(~ 1 + weight, pen_lsl_load_covar),
                              Z = cbind(Z_grm_pen_lsl_load, Z_pen_lsl_load_group),
                              RandC = c(ncol(Z_grm_pen_lsl_load), ncol(Z_pen_lsl_load_group)),
                              calc.like = TRUE)

pen_lsl_load_null <- hglm(y = pen_lsl_load$data$X6,
                             X = model.matrix(~ 1 + weight, pen_lsl_load_covar),
                             Z = Z_pen_lsl_load_group,
                             calc.like = TRUE)

lrt(pen_lsl_load_model, pen_lsl_load_null)

cage_lsl_load_model <- hglm(y = cage_lsl_load$data$X6,
                               X = model.matrix(~ 1 + weight, cage_lsl_load_covar),
                               Z = Z_grm_cage_lsl_load,
                               calc.like = TRUE)


lrt(cage_bovans_load_model)


pen_lsl_load_model$SummVC2
get_h2(pen_lsl_load_model)
get_group_ratio(pen_lsl_load_model)

cage_lsl_load_model$SummVC2
get_h2(cage_lsl_load_model)


pen_bovans_weight_model <- hglm(y = pen_bovans_weight$data$X6,
                                X = model.matrix(~ 1, pen_bovans_weight_covar),
                                Z = cbind(Z_grm_pen_bovans_weight, Z_pen_bovans_weight_group),
                                RandC = c(ncol(Z_grm_pen_bovans_weight), ncol(Z_pen_bovans_weight_group)),
                                calc.lik = TRUE)

pen_bovans_weight_null <- hglm(y = pen_bovans_weight$data$X6,
                               X = model.matrix(~ 1, pen_bovans_weight_covar),
                               Z = Z_pen_bovans_weight_group,
                               calc.like = TRUE)

lrt(pen_bovans_weight_model, pen_bovans_weight_null)

cage_bovans_weight_model <- hglm(y = cage_bovans_weight$data$X6,
                                 X = model.matrix(~ 1, cage_bovans_weight_covar),
                                 Z = Z_grm_cage_bovans_weight,
                                 calc.like = TRUE)

lrt(cage_bovans_weight_model)


pen_bovans_weight_model$SummVC2
get_h2(pen_bovans_weight_model)
get_group_ratio(pen_bovans_weight_model)

cage_bovans_weight_model$SummVC2
get_h2(cage_bovans_weight_model)


pen_lsl_weight_model <- hglm(y = pen_lsl_weight$data$X6,
                             X = model.matrix(~ 1, pen_lsl_weight_covar),
                             Z = cbind(Z_grm_pen_lsl_weight, Z_pen_lsl_weight_group),
                             RandC = c(ncol(Z_grm_pen_lsl_weight), ncol(Z_pen_lsl_weight_group)),
                             calc.lik = TRUE)

pen_lsl_weight_null <- hglm(y = pen_lsl_weight$data$X6,
                            X = model.matrix(~ 1, pen_lsl_weight_covar),
                            Z = Z_pen_lsl_weight_group,
                            calc.like = TRUE)

lrt(pen_lsl_weight_model, pen_lsl_weight_null)

cage_lsl_weight_model <- hglm(y = cage_lsl_weight$data$X6,
                              X = model.matrix(~ 1, cage_lsl_weight_covar),
                              Z = Z_grm_cage_lsl_weight,
                              calc.like = TRUE)

lrt(cage_lsl_weight_model)


pen_lsl_weight_model$SummVC2
get_h2(pen_lsl_weight_model)
get_group_ratio(pen_lsl_weight_model)

cage_lsl_weight_model$SummVC2
get_h2(cage_lsl_weight_model)



results <- data.frame(name = rep(c("bone strength", "body weight"), each = 2),
                      cross = rep(c("Bovans", "LSL"), 2),
                      h2_cage = c(get_h2(cage_bovans_load_model),
                                  get_h2(cage_lsl_load_model),
                                  get_h2(cage_bovans_weight_model),
                                  get_h2(cage_lsl_weight_model)),
                      lrt_cage = c(lrt(cage_bovans_load_model)$p.value,
                                   lrt(cage_lsl_load_model)$p.value,
                                   lrt(cage_bovans_weight_model)$p.value,
                                   lrt(cage_lsl_weight_model)$p.value),
                      h2_pen = c(get_h2(pen_bovans_load_model),
                                 get_h2(pen_lsl_load_model),
                                 get_h2(pen_bovans_weight_model),
                                 get_h2(pen_lsl_weight_model)),
                      lrt_pen = c(lrt(pen_bovans_load_model)$p.value,
                                  lrt(pen_lsl_load_model)$p.value,
                                  lrt(pen_bovans_weight_model)$p.value,
                                  lrt(pen_lsl_weight_model)$p.value))


write.csv(results,
          file = "tables/table_bone_strength_weight_h2_crosses_separate.csv",
          quote = TRUE,
          row.names = FALSE)
