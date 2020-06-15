

## Prepare files for running conditional GWAS on body weight locus on chr4

library(assertthat)
library(dplyr)
library(readr)

source("R/preparation_helper_functions.R")

gwas <- readRDS("outputs/gwas.Rds")

bw_chr4 <- filter(gwas, scan_name == "all_weight" & chr == 4)


peak_marker <- bw_chr4[which.min(bw_chr4$p_wald),]$rs



## Perform linear regression on peak marker, and save modified phenotype file


pheno <- readRDS("outputs/pheno.Rds")
geno_all <- readRDS("outputs/geno_all.Rds")

all_genotyped <- filter(pheno,
                        animal_id %in% geno_all$individual)

assert_that(identical(all_genotyped$animal_id,
                      geno_all$individual))


model <- lm(all_genotyped$weight ~
                all_genotyped$breed +
                all_genotyped$cage.pen + 
                geno_all$Gga_rs14490981)


all_genotyped$residual_weight <- residuals(model)


residual_fam <- fam(all_genotyped, "residual_weight")

write_plink(residual_fam, "gwas/all_residual_weight.fam")

