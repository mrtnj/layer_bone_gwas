

## Genomic models with hglm

library(assertthat)
library(dplyr)
library(hglm)
library(readr)

source("R/hglm_helper_functions.R")


source("R/hglm_gwas_tga_prepare_data.R")


pen_models <- mapply(function(data, covar, Z_grm, Z_group) {
    
    hglm(y = data$data$X6,
         X = model.matrix(~ 1 + weight + breed, covar),
         Z = cbind(Z_grm, Z_group),
         RandC = c(ncol(Z_grm), ncol(Z_group)),
         calc.like = TRUE,
         maxit = 100)
    
},
data_pen,
pen_covar,
Z_grm_pen,
Z_pen_group,
SIMPLIFY = FALSE)

lapply(pen_models, get_h2)

pen_null_models <- mapply(function(data, covar, Z_group) {
    
    hglm(y = data$data$X6,
         X = model.matrix(~ 1 + weight + breed, covar),
         Z = Z_group,
         calc.like = TRUE)
    
},
data_pen,
pen_covar,
Z_pen_group,
SIMPLIFY = FALSE)


lrt_pen <- mapply(function(model, null) lrt(model, null),
                  pen_models,
                  pen_null_models,
                  SIMPLIFY = FALSE)


cage_models <- mapply(function(data, covar, Z_grm, Z_group) {
    
    hglm(y = data$data$X6,
         X = model.matrix(~ 1 + weight + breed, covar),
         Z = Z_grm,
         calc.like = TRUE,
         maxit = 100)
    
},
data_cage,
cage_covar,
Z_grm_cage,
Z_cage_group,
SIMPLIFY = FALSE)




lapply(cage_models, get_h2)



lrt_cage <- lapply(cage_models, lrt)
                  

results <- data.frame(name = tga_traits,
                      h2_cage = unlist(lapply(cage_models, get_h2)),
                      lrt_cage = unlist(lapply(lrt_cage, function(x) x$p.value)),
                      h2_pen = unlist(lapply(pen_models, get_h2)),
                      lrt_pen = unlist(lapply(lrt_pen, function(x) x$p.value)),
                      stringsAsFactors = FALSE)


pretty_trait_names <- rbind(data.frame(name = c("ct_pc1", "ct_pc2", "ct_pc3"),
                                       pretty_name = c("pQCT PC1 'high density, thickness, content'",
                                                       "pQCT PC2 'long bone length'",
                                                       "pQCT PC3 'low cortical density'")),
                                       read_csv("pretty_trait_names_tga.csv"))

results <- inner_join(results, pretty_trait_names)


write.csv(results,
          file = "tables/table_bone_phenotype_h2.csv",
          quote = TRUE,
          row.names = FALSE)
