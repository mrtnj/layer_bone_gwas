

## Genomic models with hglm

library(assertthat)
library(dplyr)
library(hglm)
library(readr)



## What difference does a group variance component make to the pen model? Not much.

pen_load <- read_delim("gwas/pen_load_N/pen_load_N.fam", col_names = FALSE, delim = " ")

pen_missing <- is.na(pen_load$X6)

pen_load <- pen_load[!pen_missing,]


## Load GRM from Gemma and decompose

pen_GRM <- read_tsv("gwas/output/pen_grm.sXX.txt", col_names = FALSE)

svd_pen <- svd(pen_GRM[!pen_missing, !pen_missing])

Z_grm_pen <- svd_pen$u %*% diag(sqrt(svd_pen$d))


## Get group data

pheno <- readRDS("outputs/pheno.Rds")

pen_covar <- filter(pheno,
                    animal_id %in% pen_load$X1 &
                    !is.na(load_N))[, c("animal_id", "group", "weight")]

pen_covar <- pen_covar[order(pen_covar$animal_id),]

assert_that(identical(pen_covar$animal_id, pen_load$X1))

Z_pen_group <- model.matrix(~ factor(group), pen_covar)


## Model

model <- hglm(y = pen_load$X6,
              X = model.matrix(~ 1 + weight, pen_covar),
              Z = cbind(Z_grm_pen, Z_pen_group),
              RandC = c(ncol(Z_grm_pen), ncol(Z_pen_group)))


h2 <- model$varRanef[1] / (model$varFix + sum(model$varRanef))
ratio_group <- model$varRanef[2] / (model$varFix + sum(model$varRanef))


## Model without pen

model_without <- hglm(y = pen_load$X6,
                      X = model.matrix(~ 1 + weight, pen_covar),
                      Z = Z_grm_pen)

h2_without <- model_without$varRanef[1] / (model_without$varFix + sum(model_without$varRanef))


