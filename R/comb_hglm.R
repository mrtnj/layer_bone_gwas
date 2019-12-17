

## Genomic models with hglm

library(assertthat)
library(dplyr)
library(hglm)
library(readr)



## What difference does a group variance component make to the pen model? Not much.

pen_comb <- read_delim("gwas/pen_comb_adj/pen_comb_adj.fam", col_names = FALSE, delim = " ")

pen_missing <- is.na(pen_comb$X6)

pen_comb <- pen_comb[!pen_missing,]


## Load GRM from Gemma and decompose

pen_GRM <- read_tsv("gwas/output/pen_grm.sXX.txt", col_names = FALSE)

## svd_pen <- svd(pen_GRM[!pen_missing, !pen_missing])
##Z_grm_pen <- svd_pen$u %*% diag(sqrt(svd_pen$d))

L <- t(chol(pen_GRM[!pen_missing, !pen_missing]))

Z_grm_pen <- diag(nrow(pen_comb)) %*% L

## Get group data

pheno <- readRDS("outputs/pheno.Rds")

pen_covar <- filter(pheno,
                    animal_id %in% pen_comb$X1 &
                    !is.na(comb_g))[, c("animal_id", "group", "weight")]

pen_covar <- pen_covar[order(pen_covar$animal_id),]

assert_that(identical(pen_covar$animal_id, pen_comb$X1))

Z_pen_group <- model.matrix(~ factor(group), pen_covar)


## Incidence matrix for conspecifics

Z_pen_conspecific <- matrix(0,
                            nrow = nrow(pen_comb),
                            ncol = nrow(pen_comb))

colnames(Z_pen_conspecific) <- rownames(Z_pen_conspecific) <- pen_covar$animal_id

Z_grm_conspecific <- Z_pen_conspecific %*% L

for (ind_ix in 1:nrow(pen_covar)) {

    group_members <- pen_covar$animal_id[pen_covar$group == pen_covar$group[ind_ix]]
    conspecifics <- setdiff(group_members, pen_covar$animal_id[ind_ix])

    Z_pen_conspecific[ind_ix, which(colnames(Z_pen_conspecific) %in% conspecifics)] <- 1
    Z_pen_conspecific[which(colnames(Z_pen_conspecific) %in% conspecifics), ind_ix] <- 1
}


## Model

model <- hglm(y = pen_comb$X6,
              X = model.matrix(~ 1 + weight, pen_covar),
              Z = cbind(Z_grm_pen, Z_pen_group, Z_pen_conspecific),
              RandC = c(ncol(Z_grm_pen), ncol(Z_pen_group), ncol(Z_grm_conspecific)))


h2 <- model$varRanef[1] / (model$varFix + sum(model$varRanef))
ratio_group <- model$varRanef[2] / (model$varFix + sum(model$varRanef))
ratio_conspecifics <- model$varRanef[3] / (model$varFix + sum(model$varRanef))


