## Get data

pen_load <- read_data("gwas/pen_load_N/pen_load_N.fam")
cage_load <- read_data("gwas/cage_load_N/cage_load_N.fam")
all_load <- read_data("gwas/all_load_N/all_load_N.fam")

pen_weight <- read_data("gwas/pen_weight/pen_weight.fam")
cage_weight <- read_data("gwas/cage_weight/cage_weight.fam")
all_weight <- read_data("gwas/all_weight/all_weight.fam")


## Load GRM from Gemma and decompose

pen_GRM <- read_tsv("gwas/output/pen_grm.sXX.txt", col_names = FALSE)
cage_GRM <- read_tsv("gwas/output/cage_grm.sXX.txt", col_names = FALSE)
all_GRM <- read_tsv("gwas/output/all_grm.sXX.txt", col_names = FALSE)

Z_grm_pen_load <- decompose_grm(pen_GRM,
                                pen_load$missing)

Z_grm_cage_load <- decompose_grm(cage_GRM,
                                 cage_load$missing)

Z_grm_all_load <- decompose_grm(all_GRM,
                                all_load$missing)


Z_grm_pen_weight <- decompose_grm(pen_GRM,
                                  pen_weight$missing)

Z_grm_cage_weight <- decompose_grm(cage_GRM,
                                   cage_weight$missing)

Z_grm_all_weight <- decompose_grm(all_GRM,
                                  all_weight$missing)




## Group covariate and incidence matrix

pheno <- readRDS("outputs/pheno.Rds")

pen_load_covar <- get_covar(pheno,
                            pen_load$data,
                            "load_N")

cage_load_covar <- get_covar(pheno,
                             cage_load$data,
                             "load_N")

all_load_covar <- get_covar(pheno,
                            all_load$data,
                            "load_N")


pen_weight_covar <- get_covar(pheno,
                              pen_weight$data,
                              "weight")

cage_weight_covar <- get_covar(pheno,
                               cage_weight$data,
                               "weight")

all_weight_covar <- get_covar(pheno,
                              all_weight$data,
                              "weight")

all_load_covar$group_cages_combined <- ifelse(all_load_covar$group > 2000,
                                              2000,
                                              all_load_covar$group)

all_weight_covar$group_cages_combined <- ifelse(all_weight_covar$group > 2000,
                                                2000,
                                                all_weight_covar$group)

Z_pen_load_group <- model.matrix(~ factor(group), pen_load_covar)

Z_cage_load_group <- model.matrix(~ factor(group), cage_load_covar)

Z_all_load_group <- model.matrix(~ factor(group_cages_combined), all_load_covar)

Z_pen_weight_group <- model.matrix(~ factor(group), pen_weight_covar)

Z_cage_weight_group <- model.matrix(~ factor(group), cage_weight_covar)

Z_all_weight_group <- model.matrix(~ factor(group_cages_combined), all_weight_covar)



## Genotypes

geno_pen <- read_delim("gwas/pen.raw",
                       delim = " ")

geno_pen_load <- geno_pen[match(pen_load$data$X1, geno_pen$FID),]
assert_that(identical(geno_pen_load$FID, pen_load$data$X1))

geno_pen_weight <- geno_pen[match(pen_weight$data$X1, geno_pen$FID),]
assert_that(identical(geno_pen_weight$FID, pen_weight$data$X1))


geno_cage <- read_delim("gwas/cage.raw",
                        delim = " ")

geno_cage_load <- geno_cage[match(cage_load$data$X1, geno_cage$FID),]
assert_that(identical(geno_cage_load$FID, cage_load$data$X1))

geno_cage_weight <- geno_cage[match(cage_weight$data$X1, geno_cage$FID),]
assert_that(identical(geno_cage_weight$FID, cage_weight$data$X1))


geno_all <- read_delim("gwas/all.raw",
                       delim = " ")

geno_all_load <- geno_all[match(all_load$data$X1, geno_all$FID),]
assert_that(identical(geno_all_load$FID, all_load$data$X1))

geno_all_weight <- geno_all[match(all_weight$data$X1, geno_all$FID),]
assert_that(identical(geno_all_weight$FID, all_weight$data$X1))


snps_pruned_pen_load <- prune_snp_matrix(geno_pen_load[, -(1:6)])
snps_pruned_pen_weight <- prune_snp_matrix(geno_pen_weight[, -(1:6)])

snps_pruned_cage_load <- prune_snp_matrix(geno_cage_load[, -(1:6)])
snps_pruned_cage_weight <- prune_snp_matrix(geno_cage_weight[, -(1:6)])

snps_pruned_all_load <- prune_snp_matrix(geno_all_load[, -(1:6)])
snps_pruned_all_weight <- prune_snp_matrix(geno_all_weight[, -(1:6)])
