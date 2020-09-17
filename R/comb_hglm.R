

## Genomic models with hglm

library(assertthat)
library(dplyr)
library(hglm)
library(readr)


source("R/hglm_helper_functions.R")



## What difference does a group variance component make to the pen model? Not much.

pen_comb <- read_data("gwas/pen_comb_g/pen_comb_g.fam")

cage_comb <- read_data("gwas/cage_comb_g/cage_comb_g.fam")

# pen_missing <- is.na(pen_comb$data$X6)
# 
# pen_comb <- pen_comb[!pen_missing,]


## Load GRM from Gemma and decompose

pen_GRM <- read_tsv("gwas/output/pen_grm.sXX.txt", col_names = FALSE)

cage_GRM <- read_tsv("gwas/output/cage_grm.sXX.txt", col_names = FALSE)

Z_grm_pen <- decompose_grm(pen_GRM,
                           pen_comb$missing)

Z_grm_cage <- decompose_grm(cage_GRM,
                            cage_comb$missing)

## svd_pen <- svd(pen_GRM[!pen_missing, !pen_missing])
##Z_grm_pen <- svd_pen$u %*% diag(sqrt(svd_pen$d))

L_pen <- t(chol(pen_GRM[!pen_comb$missing, !pen_comb$missing]))

L_cage <- t(chol(cage_GRM[!cage_comb$missing, !cage_comb$missing]))
# 
# Z_grm_pen <- diag(nrow(pen_comb)) %*% L

## Get group data

pheno <- readRDS("outputs/pheno.Rds")

pen_covar <- get_covar(pheno,
                       pen_comb$data,
                       "comb_g")

cage_covar <- get_covar(pheno,
                        cage_comb$data,
                        "comb_g")

# pen_covar <- filter(pheno,
#                     animal_id %in% pen_comb$X1 &
#                     !is.na(comb_g))[, c("animal_id", "group", "weight")]
# 
# pen_covar <- pen_covar[order(pen_covar$animal_id),]
# 
# assert_that(identical(pen_covar$animal_id, pen_comb$X1))

Z_pen_group <- model.matrix(~ factor(group), pen_covar)

Z_cage_group <- model.matrix(~ factor(group), cage_covar)


## Model

model_pen <- hglm(y = pen_comb$data$X6,
                  X = model.matrix(~ 1 + weight + breed, pen_covar),
                  Z = cbind(Z_grm_pen, Z_pen_group),
                  RandC = c(ncol(Z_grm_pen), ncol(Z_pen_group)),
                  conv = 1e-8)


h2 <- get_h2(model_pen)
ratio_group <- get_group_ratio(model_pen)



model_cage <- hglm(y = cage_comb$data$X6,
                   X = model.matrix(~ 1 + weight + breed, cage_covar),
                   Z = cbind(Z_grm_cage, Z_cage_group),
                   RandC = c(ncol(Z_grm_cage), ncol(Z_cage_group)),
                   conv = 1e-8)


h2_cage <- get_h2(model_cage)
ratio_group_cage <- get_group_ratio(model_cage)








## Social model that includes genetic effects of conspecifics

## Incidence matrix for conspecifics

get_conspecific_matrix <- function(data,
                                   covar,
                                   L) {

    Z_conspecific <- matrix(0,
                            nrow = nrow(data$data),
                            ncol = nrow(data$data))
    
    colnames(Z_conspecific) <- rownames(Z_conspecific) <- covar$animal_id
    
    for (ind_ix in 1:nrow(covar)) {
        
        group_members <- covar$animal_id[covar$group == covar$group[ind_ix]]
        conspecifics <- setdiff(group_members, covar$animal_id[ind_ix])
        
        Z_conspecific[ind_ix, which(colnames(Z_conspecific) %in% conspecifics)] <- 1
        Z_conspecific[which(colnames(Z_conspecific) %in% conspecifics), ind_ix] <- 1
    }
    
    Z_grm_conspecific <- Z_conspecific %*% L
    
    Z_grm_conspecific
}

Z_pen_conspecific <- get_conspecific_matrix(pen_comb,
                                            pen_covar,
                                            L_pen)

Z_cage_conspecific <- get_conspecific_matrix(cage_comb,
                                             cage_covar,
                                             L_cage)

## Model

model_social <- hglm(y = pen_comb$data$X6,
                     X = model.matrix(~ 1 + weight, pen_covar),
                     Z = cbind(Z_grm_pen, Z_pen_group, Z_grm_conspecific),
                     RandC = c(ncol(Z_grm_pen), ncol(Z_pen_group), ncol(Z_grm_conspecific)),
                     conv = 1e-8)


h2 <- get_var_ratio(model_social, 1)
ratio_group <- get_var_ratio(model_social, 2)
ratio_conspecifics <- get_var_ratio(model_social, 3)



model_social_cage <- hglm(y = cage_comb$data$X6,
                     X = model.matrix(~ 1 + weight, cage_covar),
                     Z = cbind(Z_grm_cage, Z_cage_group, Z_cage_conspecific),
                     RandC = c(ncol(Z_grm_cage), ncol(Z_cage_group), ncol(Z_cage_conspecific)),
                     conv = 1e-8)


h2_cage_social <- get_var_ratio(model_social_cage, 1)
ratio_group_cage_social <- get_var_ratio(model_social_cage, 2)
ratio_conspecifics_cage_social <- get_var_ratio(model_social_cage, 3)
