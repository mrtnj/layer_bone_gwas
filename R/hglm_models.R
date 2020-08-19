

## Genomic models with hglm

library(assertthat)
library(dplyr)
library(hglm)
library(readr)

source("R/hglm_helper_functions.R")


source("R/hglm_gwas_prepare_data.R")



pen_load <- hglm(y = pen_load$data$X6,
                 X = model.matrix(~ 1 + weight + breed, pen_load_covar),
                 Z = cbind(Z_grm_pen_load, Z_pen_load_group),
                 RandC = c(ncol(Z_grm_pen_load), ncol(Z_pen_load_group)))

cage_load <- hglm(y = cage_load$data$X6,
                  X = model.matrix(~ 1 + weight + breed, cage_load_covar),
                  Z = Z_grm_cage_load,
                  RandC = ncol(Z_grm_cage_load))

pen_load$SummVC2
get_h2(pen_load)
get_group_ratio(pen_load)

cage_load$SummVC2
get_h2(cage_load)



pen_weight <- hglm(y = pen_weight$data$X6,
                   X = model.matrix(~ 1 + breed, pen_weight_covar),
                   Z = cbind(Z_grm_pen_weight, Z_pen_weight_group),
                   RandC = c(ncol(Z_grm_pen_weight), ncol(Z_pen_weight_group)))

cage_weight <- hglm(y = cage_weight$data$X6,
                    X = model.matrix(~ 1 + breed, cage_weight_covar),
                    Z = Z_grm_cage_weight,
                    RandC = ncol(Z_grm_cage_weight))


pen_weight$SummVC2
get_h2(pen_weight)
get_group_ratio(pen_weight)

cage_weight$SummVC2
get_h2(cage_weight)



