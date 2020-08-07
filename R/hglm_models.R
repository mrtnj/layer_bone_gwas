

## Genomic models with hglm

library(assertthat)
library(dplyr)
library(hglm)
library(readr)



## What difference does a group variance component make to the pen model? Not much.

pen_load <- read_delim("gwas/pen_load_N/pen_load_N.fam", col_names = FALSE, delim = " ")

pen_missing <- is.na(pen_load$X6)

pen_load <- pen_load[!pen_missing,]


all_load <- read_delim("gwas/all_load_N/all_load_N.fam", col_names = FALSE, delim = " ")

all_missing <- is.na(all_load$X6)

all_load <- all_load[!all_missing,]


## Load GRM from Gemma and decompose

pen_GRM <- read_tsv("gwas/output/pen_grm.sXX.txt", col_names = FALSE)

all_GRM <- read_tsv("gwas/output/all_grm.sXX.txt", col_names = FALSE)

decompose_grm <- function(grm, missing) {
    svd <- svd(grm[!missing, !missing])
    svd$u %*% diag(sqrt(svd$d))
}

Z_grm_pen <- decompose_grm(pen_GRM,
                           pen_missing)

Z_grm_all <- decompose_grm(all_GRM,
                           all_missing)


## Get group data

pheno <- readRDS("outputs/pheno.Rds")

get_covar <- function(pheno,
                      load_data) {

    covar <- filter(pheno,
                        animal_id %in% load_data$X1 &
                            !is.na(load_N))[, c("animal_id", "group", "weight")]
    
    covar <- covar[order(covar$animal_id),]
    
    assert_that(identical(covar$animal_id, load_data$X1))

    covar
}

pen_covar <- get_covar(pheno,
                       pen_load)

all_covar <- get_covar(pheno,
                       all_load)

all_covar$group_cages_combined <- ifelse(all_covar$group > 2000,
                                         2000,
                                         all_covar$group)

Z_pen_group <- model.matrix(~ factor(group), pen_covar)

Z_all_group <- model.matrix(~ factor(group), all_covar)

Z_all_group_cages_combined <- model.matrix(~ factor(group_cages_combined), all_covar)


get_h2 <- function(model) {
    model$varRanef[1] / (model$varFix + sum(model$varRanef))
}

get_group_ratio <- function(model) {
    model$varRanef[2] / (model$varFix + sum(model$varRanef))
}

## Model of pen animals with random pen group

model_pen_random <- hglm(y = pen_load$X6,
                         X = model.matrix(~ 1 + weight, pen_covar),
                         Z = cbind(Z_grm_pen, Z_pen_group),
                         RandC = c(ncol(Z_grm_pen), ncol(Z_pen_group)))


h2_pen_random <- get_h2(model_pen_random)
ratio_group <- get_group_ratio(model_pen_random)


## Model of pen animals without random pen group

model_pen_without <- hglm(y = pen_load$X6,
                          X = model.matrix(~ 1 + weight, pen_covar),
                          Z = Z_grm_pen)

h2_pen_without <- get_h2(model_pen_without)


## All animals with random pen/cage group

model_all_random <- hglm(y = all_load$X6,
                         X = model.matrix(~ 1 + weight, all_covar),
                         Z = cbind(Z_grm_all, Z_all_group),
                         RandC = c(ncol(Z_grm_all), ncol(Z_all_group)))


h2_all_random <- get_h2(model_all_random)
ratio_group_all_random <- get_group_ratio(model_all_random)


## All animals with random pen group, cages combined into one group;

model_all_random_cages_combined <- hglm(y = all_load$X6,
                                        X = model.matrix(~ 1 + weight, all_covar),
                                        Z = cbind(Z_grm_all, Z_all_group_cages_combined),
                                        RandC = c(ncol(Z_grm_all), ncol(Z_all_group_cages_combined)))


h2_all_random_cages_combined <- get_h2(model_all_random_cages_combined)
ratio_group_all_random_cages_combined <- get_group_ratio(model_all_random_cages_combined)


## All animals without random pen/cage group

model_all_without <- hglm(y = all_load$X6,
                          X = model.matrix(~ 1 + weight, all_covar),
                          Z = Z_grm_all)

h2_all_without <- get_h2(model_all_without)
