

## Perform GWAS with hglm and extra variance components

library(assertthat)
library(dplyr)
library(hglm)
library(readr)


source("R/hglm_helper_functions.R")


pen_load <- read_delim("gwas/pen_load_N/pen_load_N.fam", col_names = FALSE, delim = " ")

pen_missing <- is.na(pen_load$X6)

pen_load <- pen_load[!pen_missing,]


all_load <- read_delim("gwas/all_load_N/all_load_N.fam", col_names = FALSE, delim = " ")

all_missing <- is.na(all_load$X6)

all_load <- all_load[!all_missing,]


## Load GRM from Gemma and decompose

pen_GRM <- read_tsv("gwas/output/pen_grm.sXX.txt", col_names = FALSE)

all_GRM <- read_tsv("gwas/output/all_grm.sXX.txt", col_names = FALSE)

Z_grm_pen <- decompose_grm(pen_GRM,
                           pen_missing)

Z_grm_all <- decompose_grm(all_GRM,
                           all_missing)


## Get group data

pheno <- readRDS("outputs/pheno.Rds")

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



## Genotypes

geno_pen <- read_delim("gwas/pen.raw",
                       delim = " ")

geno_pen <- geno_pen[!pen_missing,]

assert_that(identical(geno_pen$FID, pen_load$X1))

geno_all <- read_delim("gwas/all.raw",
                       delim = " ")

geno_all <- geno_all[!all_missing,]

assert_that(identical(geno_all$FID, all_load$X1))

snps_pruned_pen <- prune_snp_matrix(geno_pen)
snps_pruned_all <- prune_snp_matrix(geno_all)


## Run GWAS

run_gwas <- function(pheno,
                     X,
                     Z,
                     RandC,
                     snp_matrix) {
    
    n_ind <- length(pheno)
    n_snp <- ncol(snp_matrix)

    ## Fit baseline model

    model_baseline <- hglm(y = pheno,
                         X = X,
                         Z = Z,
                         RandC = RandC,
                         X.disp = X.disp)

    
    ratio <- model_baseline$varRanef/model_pen_random$varFix

    V <- RepeatABEL::constructV(Z,
                                RandC,
                                ratio)

    eigV <- eigen(V)

    transformation_matrix <- diag(1/sqrt(eigV$values)) %*% t(eigV$vectors)
    
    transformed_y <- transformation_matrix %*% pheno
    transformed_X <- transformation_matrix %*% X


    ## Null model
    
    qr0 <- qr(transformed_X)
    
    est0 <- qr.coef(qr0, transformed_y)

    null_residual <- transformed_y - transformed_X %*% est0

    RSS_null <- sum(null_residual^2)/n_ind


    ## SNP models
    
    estimates <- numeric(n_snp)
    LRT <- numeric(n_snp)
    p <- numeric(n_snp)

    for (snp_ix in 1:n_snp) {
        
        transformed_snp <- transformation_matrix %*% as.matrix(snp_matrix[, snp_ix])
        
        X1 <- cbind(transformed_snp, transformed_X)
        qr1 <- qr(X1)
        est1 <- qr.coef(qr1, transformed_y)
        residual1 <- transformed_y - X1 %*% est1
        RSS1 <- sum(residual1^2)/n_ind
        
        estimates[snp_ix] <- est1[1]
        LRT[snp_ix] <- -n_ind * (log(RSS1) - log(RSS_null))
        p[snp_ix] <- 1 - pchisq(LRT[snp_ix],
                                df = 1)
    }
    
    data.frame(marker_id = colnames(snp_matrix),
               estimates,
               LRT,
               p)
}

gwas_pen_random <- run_gwas(pen_load$X6,
                            model.matrix(~ 1 + weight, pen_covar),
                            cbind(Z_grm_pen, Z_pen_group),
                            c(ncol(Z_grm_pen), ncol(Z_pen_group)),
                            snps_pruned_pen)



gwas_all_random <- run_gwas(all_load$X6,
                            model.matrix(~ 1 + weight + cage.pen, all_covar),
                            cbind(Z_grm_all, Z_all_group),
                            c(ncol(Z_grm_all), ncol(Z_all_group)),
                            snps_pruned_all)

