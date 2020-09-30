

## Genomic models with hglm

library(assertthat)
library(dplyr)
library(hglm)
library(readr)

source("R/hglm_helper_functions.R")


source("R/hglm_gwas_prepare_data.R")



pen_load_model <- hglm(y = pen_load$data$X6,
                 X = model.matrix(~ 1 + weight + breed, pen_load_covar),
                 Z = cbind(Z_grm_pen_load, Z_pen_load_group),
                 RandC = c(ncol(Z_grm_pen_load), ncol(Z_pen_load_group)),
                 calc.like = TRUE)

pen_load_null <- hglm(y = pen_load$data$X6,
                      X = model.matrix(~ 1 + weight + breed, pen_load_covar),
                      Z = Z_pen_load_group,
                      calc.like = TRUE)

lrt(pen_load_model, pen_load_null)

cage_load_model <- hglm(y = cage_load$data$X6,
                  X = model.matrix(~ 1 + weight + breed, cage_load_covar),
                  Z = Z_grm_cage_load,
                  calc.like = TRUE)


lrt(cage_load_model)

pen_load_model$SummVC2
get_h2(pen_load_model)
get_group_ratio(pen_load_model)

cage_load_model$SummVC2
get_h2(cage_load_model)



pen_weight_model <- hglm(y = pen_weight$data$X6,
                         X = model.matrix(~ 1 + breed, pen_weight_covar),
                         Z = cbind(Z_grm_pen_weight, Z_pen_weight_group),
                         RandC = c(ncol(Z_grm_pen_weight), ncol(Z_pen_weight_group)),
                         calc.lik = TRUE)

pen_weight_null <- hglm(y = pen_weight$data$X6,
                        X = model.matrix(~ 1 + breed, pen_weight_covar),
                        Z = Z_pen_weight_group,
                        calc.like = TRUE)

lrt(pen_weight_model, pen_weight_null)

cage_weight_model <- hglm(y = cage_weight$data$X6,
                          X = model.matrix(~ 1 + breed, cage_weight_covar),
                          Z = Z_grm_cage_weight,
                          calc.like = TRUE)

lrt(cage_weight_model)


pen_weight_model$SummVC2
get_h2(pen_weight_model)
get_group_ratio(pen_weight_model)

cage_weight_model$SummVC2
get_h2(cage_weight_model)



results <- data.frame(name = c("bone strength", "body weight"),
                      h2_cage = c(get_h2(cage_load_model), get_h2(cage_weight_model)),
                      lrt_cage = c(lrt(cage_load_model)$p.value, lrt(cage_weight_model)$p.value),
                      h2_pen = c(get_h2(pen_load_model), get_h2(pen_weight_model)),
                      lrt_pen = c(lrt(pen_load_model)$p.value, lrt(pen_weight_model)$p.value))


write.csv(results,
          file = "tables/table_bone_strength_weight_h2.csv",
          quote = TRUE,
          row.names = FALSE)
