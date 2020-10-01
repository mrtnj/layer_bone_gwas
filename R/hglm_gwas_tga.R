

## Perform GWAS on selected pQCT/TGA traits

library(assertthat)
library(dplyr)
library(hglm)
library(readr)


source("R/hglm_helper_functions.R")


source("R/hglm_gwas_tga_prepare_data.R")


h2_results <- read_csv("tables/table_bone_phenotype_h2.csv")

cage_traits <- filter(h2_results, lrt_cage < 0.05)$name
pen_traits <- filter(h2_results, lrt_pen < 0.05)$name
all_traits <- intersect(cage_traits, pen_traits)

cage_trait_ix <- match(cage_traits, tga_traits)
pen_trait_ix <- match(pen_traits, tga_traits)
all_trait_ix <- match(all_traits, tga_traits)

gwas_pen <- mapply(function(data, covar, Z_grm, Z_group, snps, trait_name) {
    gwas <- run_gwas(data$data$X6,
                     model.matrix(~ 1 + weight + breed, covar),
                     cbind(Z_grm, Z_group),
                     c(ncol(Z_grm), ncol(Z_group)),
                     snps)
    gwas$trait <- trait_name
    gwas
},
data_pen[pen_trait_ix],
pen_covar[pen_trait_ix],
Z_grm_pen[pen_trait_ix],
Z_pen_group[pen_trait_ix],
snps_pruned_pen[pen_trait_ix],
pen_traits,
SIMPLIFY = FALSE)
    
gwas_pen_df <- Reduce(rbind, gwas_pen)


saveRDS(gwas_pen_df,
        file = "gwas/hglm_gwas_pen_bone_phenotypes.Rds")


gwas_cage <- mapply(function(data, covar, Z_grm, snps, trait_name) {
    gwas <- run_gwas(data$data$X6,
                     model.matrix(~ 1 + weight + breed, covar),
                     Z_grm,
                     ncol(Z_grm),
                     snps)
    gwas$trait <- trait_name
    gwas
},
data_cage[cage_trait_ix],
cage_covar[cage_trait_ix],
Z_grm_cage[cage_trait_ix],
snps_pruned_cage[cage_trait_ix],
cage_traits,
SIMPLIFY = FALSE)


gwas_cage_df <- Reduce(rbind, gwas_cage)


saveRDS(gwas_cage_df,
        file = "gwas/hglm_gwas_cage_bone_phenotypes.Rds")



gwas_all <- mapply(function(data, covar, Z_grm, Z_group, snps, trait_name) {
    gwas <- run_gwas(data$data$X6,
                     model.matrix(~ 1 + weight + breed + cage.pen, covar),
                     cbind(Z_grm, Z_group),
                     c(ncol(Z_grm), ncol(Z_group)),
                     snps)
    gwas$trait <- trait_name
    gwas
},
data_all[all_trait_ix],
all_covar[all_trait_ix],
Z_grm_all[all_trait_ix],
Z_all_group[all_trait_ix],
snps_pruned_all[all_trait_ix],
all_traits,
SIMPLIFY = FALSE)


gwas_all_df <- Reduce(rbind, gwas_all)

saveRDS(gwas_all,
        file = "gwas/hglm_gwas_all_bone_phenotypes.Rds")



