tga_traits <- c("ct_pc1", "ct_pc2", "ct_pc3",
                "CO2Lost_CB", "CO2Lost_MB",
                "CO3_over_Phosphates_CB", "CO3_over_Phosphates_MB",
                "Mineral_CB", "Mineral_MB",
                "OMLost_CB", "OMLost_MB",    
                "Phosphates_CB", "Phosphates_MB",
                "Phosphates_over_OM_CB", "Phosphates_over_OM_MB",
                "WaterLost_CB", "WaterLost_MB")



## Get data

data_pen <- lapply(tga_traits,
                   function(trait)
                       read_data(paste("gwas/pen_",
                                       trait,
                                       "/pen_",
                                       trait,
                                       ".fam",
                                       sep = "")))

data_cage <- lapply(tga_traits,
                    function(trait)
                        read_data(paste("gwas/cage_",
                                        trait,
                                        "/cage_",
                                        trait,
                                        ".fam",
                                        sep = "")))

data_all <- lapply(tga_traits,
                   function(trait)
                       read_data(paste("gwas/all_",
                                       trait,
                                       "/all_",
                                       trait,
                                       ".fam",
                                       sep = "")))



## Load GRM from Gemma and decompose

pen_GRM <- read_tsv("gwas/output/pen_grm.sXX.txt", col_names = FALSE)
cage_GRM <- read_tsv("gwas/output/cage_grm.sXX.txt", col_names = FALSE)
all_GRM <- read_tsv("gwas/output/all_grm.sXX.txt", col_names = FALSE)


Z_grm_pen <- lapply(data_pen,
                    function(data)
                        decompose_grm(pen_GRM,
                                      data$missing))

Z_grm_cage <- lapply(data_cage,
                    function(data)
                        decompose_grm(cage_GRM,
                                      data$missing))

Z_grm_all <- lapply(data_all,
                    function(data)
                        decompose_grm(all_GRM,
                                      data$missing))



## Group covariate and incidence matrix

pheno <- readRDS("outputs/pheno.Rds")


pen_covar <- mapply(function(data, trait) get_covar(pheno,
                                                    data$data,
                                                    trait),
                    data_pen,
                    tga_traits,
                    SIMPLIFY = FALSE)

cage_covar <- mapply(function(data, trait) get_covar(pheno,
                                                    data$data,
                                                    trait),
                    data_cage,
                    tga_traits,
                    SIMPLIFY = FALSE)

all_covar <- mapply(function(data, trait) get_covar(pheno,
                                                    data$data,
                                                    trait),
                    data_all,
                    tga_traits,
                    SIMPLIFY = FALSE)


all_covar <- lapply(all_covar,
                    function(covar) {
                        covar$group_cages_combined <- ifelse(covar$group > 2000,
                                                             2000,
                                                             covar$group)
                        covar
                    })


Z_pen_group <- lapply(pen_covar,
                      function(data) model.matrix(~ factor(group), data))

Z_cage_group <- lapply(cage_covar,
                      function(data) model.matrix(~ factor(group), data))

Z_all_group <- lapply(all_covar,
                      function(data) model.matrix(~ factor(group_cages_combined), data))




## Genotypes

geno_pen <- read_delim("gwas/pen.raw",
                       delim = " ")

geno_pen_trait <- lapply(data_pen,
                         function(data) {
                             g <- geno_pen[match(data$data$X1, geno_pen$FID),]
                             assert_that(identical(g$FID, data$data$X1))
                             g
                         })


geno_cage <- read_delim("gwas/cage.raw",
                        delim = " ")

geno_cage_trait <- lapply(data_cage,
                          function(data) {
                              g <- geno_cage[match(data$data$X1, geno_cage$FID),]
                              assert_that(identical(g$FID, data$data$X1))
                              g
                          })



geno_all <- read_delim("gwas/all.raw",
                       delim = " ")

geno_all_trait <- lapply(data_all,
                          function(data) {
                              g <- geno_all[match(data$data$X1, geno_all$FID),]
                              assert_that(identical(g$FID, data$data$X1))
                              g
                          })


snps_pruned_pen <- lapply(geno_pen_trait,
                          function(g) prune_snp_matrix(g[, -(1:6)]))

snps_pruned_cage <- lapply(geno_cage_trait,
                           function(g) prune_snp_matrix(g[, -(1:6)]))

snps_pruned_all <- lapply(geno_all_trait,
                          function(g) prune_snp_matrix(g[, -(1:6)]))

