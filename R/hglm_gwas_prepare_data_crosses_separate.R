## Prepare data for cross-separated GWAS

## Get data

pen_bovans_load <- read_data("gwas/fam_pen_bovans_load_N.fam")
cage_bovans_load <- read_data("gwas/fam_cage_bovans_load_N.fam")
all_bovans_load <- read_data("gwas/fam_all_bovans_load_N.fam")

pen_lsl_load <- read_data("gwas/fam_pen_lsl_load_N.fam")
cage_lsl_load <- read_data("gwas/fam_cage_lsl_load_N.fam")
all_lsl_load <- read_data("gwas/fam_all_lsl_load_N.fam")

pen_bovans_weight <- read_data("gwas/fam_pen_bovans_load_N.fam")
cage_bovans_weight <- read_data("gwas/fam_cage_bovans_load_N.fam")
all_bovans_weight <- read_data("gwas/fam_all_bovans_load_N.fam")

pen_lsl_weight <- read_data("gwas/fam_pen_lsl_weight.fam")
cage_lsl_weight <- read_data("gwas/fam_cage_lsl_weight.fam")
all_lsl_weight <- read_data("gwas/fam_all_lsl_weight.fam")



## Load GRM from Gemma and decompose

pen_bovans_GRM <- read_tsv("gwas/output/bovans_pen_grm.sXX.txt", col_names = FALSE)
cage_bovans_GRM <- read_tsv("gwas/output/bovans_cage_grm.sXX.txt", col_names = FALSE)
all_bovans_GRM <- read_tsv("gwas/output/all_bovans_grm.sXX.txt", col_names = FALSE)

pen_lsl_GRM <- read_tsv("gwas/output/lsl_pen_grm.sXX.txt", col_names = FALSE)
cage_lsl_GRM <- read_tsv("gwas/output/lsl_cage_grm.sXX.txt", col_names = FALSE)
all_lsl_GRM <- read_tsv("gwas/output/all_lsl_grm.sXX.txt", col_names = FALSE)


Z_grm_pen_bovans_load <- decompose_grm(pen_bovans_GRM,
                                       pen_bovans_load$missing)

Z_grm_cage_bovans_load <- decompose_grm(cage_bovans_GRM,
                                        cage_bovans_load$missing)

Z_grm_all_bovans_load <- decompose_grm(all_bovans_GRM,
                                       all_bovans_load$missing)


Z_grm_pen_lsl_load <- decompose_grm(pen_lsl_GRM,
                                    pen_lsl_load$missing)

Z_grm_cage_lsl_load <- decompose_grm(cage_lsl_GRM,
                                     cage_lsl_load$missing)

Z_grm_all_lsl_load <- decompose_grm(all_lsl_GRM,
                                    all_lsl_load$missing)



Z_grm_pen_bovans_weight <- decompose_grm(pen_bovans_GRM,
                                         pen_bovans_weight$missing)

Z_grm_cage_bovans_weight <- decompose_grm(cage_bovans_GRM,
                                          cage_bovans_weight$missing)

Z_grm_all_bovans_weight <- decompose_grm(all_bovans_GRM,
                                         all_bovans_weight$missing)


Z_grm_pen_lsl_weight <- decompose_grm(pen_lsl_GRM,
                                      pen_lsl_weight$missing)

Z_grm_cage_lsl_weight <- decompose_grm(cage_lsl_GRM,
                                       cage_lsl_weight$missing)

Z_grm_all_lsl_weight <- decompose_grm(all_lsl_GRM,
                                      all_lsl_weight$missing)



## Group covariate and incidence matrix

pheno <- readRDS("outputs/pheno.Rds")

pen_bovans_load_covar <- get_covar(pheno,
                                   pen_bovans_load$data,
                                   "load_N")

cage_bovans_load_covar <- get_covar(pheno,
                                    cage_bovans_load$data,
                                    "load_N")

all_bovans_load_covar <- get_covar(pheno,
                                   all_bovans_load$data,
                                   "load_N")


pen_lsl_load_covar <- get_covar(pheno,
                                pen_lsl_load$data,
                                "load_N")

cage_lsl_load_covar <- get_covar(pheno,
                                 cage_lsl_load$data,
                                 "load_N")

all_lsl_load_covar <- get_covar(pheno,
                                all_lsl_load$data,
                                "load_N")



pen_bovans_weight_covar <- get_covar(pheno,
                                     pen_bovans_weight$data,
                                     "load_N")

cage_bovans_weight_covar <- get_covar(pheno,
                                      cage_bovans_weight$data,
                                      "load_N")

all_bovans_weight_covar <- get_covar(pheno,
                                     all_bovans_weight$data,
                                     "load_N")


pen_lsl_weight_covar <- get_covar(pheno,
                                  pen_lsl_weight$data,
                                  "load_N")

cage_lsl_weight_covar <- get_covar(pheno,
                                   cage_lsl_weight$data,
                                   "load_N")

all_lsl_weight_covar <- get_covar(pheno,
                                  all_lsl_weight$data,
                                  "load_N")



all_bovans_load_covar$group_cages_combined <- ifelse(all_bovans_load_covar$group > 2000,
                                                     2000,
                                                     all_bovans_load_covar$group)

all_lsl_load_covar$group_cages_combined <- ifelse(all_lsl_load_covar$group > 2000,
                                                  2000,
                                                  all_lsl_load_covar$group)


all_bovans_weight_covar$group_cages_combined <- ifelse(all_bovans_weight_covar$group > 2000,
                                                     2000,
                                                     all_bovans_weight_covar$group)

all_lsl_weight_covar$group_cages_combined <- ifelse(all_lsl_weight_covar$group > 2000,
                                                    2000,
                                                    all_lsl_weight_covar$group)


Z_pen_bovans_load_group <- model.matrix(~ factor(group), pen_bovans_load_covar)

Z_cage_bovans_load_group <- model.matrix(~ factor(group), cage_bovans_load_covar)

Z_all_bovans_load_group <- model.matrix(~ factor(group_cages_combined), all_bovans_load_covar)


Z_pen_lsl_load_group <- model.matrix(~ factor(group), pen_lsl_load_covar)

Z_cage_lsl_load_group <- model.matrix(~ factor(group), cage_lsl_load_covar)

Z_all_lsl_load_group <- model.matrix(~ factor(group_cages_combined), all_lsl_load_covar)



Z_pen_bovans_weight_group <- model.matrix(~ factor(group), pen_bovans_weight_covar)

Z_cage_bovans_weight_group <- model.matrix(~ factor(group), cage_bovans_weight_covar)

Z_all_bovans_weight_group <- model.matrix(~ factor(group_cages_combined), all_bovans_weight_covar)


Z_pen_lsl_weight_group <- model.matrix(~ factor(group), pen_lsl_weight_covar)

Z_cage_lsl_weight_group <- model.matrix(~ factor(group), cage_lsl_weight_covar)

Z_all_lsl_weight_group <- model.matrix(~ factor(group_cages_combined), all_lsl_weight_covar)




## Genotypes

geno_pen_bovans <- read_delim("gwas/pen_bovans.raw",
                              delim = " ")

geno_pen_bovans_load <- geno_pen_bovans[match(pen_bovans_load$data$X1, geno_pen_bovans$FID),]
assert_that(identical(geno_pen_bovans_load$FID, pen_bovans_load$data$X1))

geno_pen_bovans_weight <- geno_pen_bovans[match(pen_bovans_weight$data$X1, geno_pen_bovans$FID),]
assert_that(identical(geno_pen_bovans_weight$FID, pen_bovans_weight$data$X1))


geno_cage_bovans <- read_delim("gwas/cage_bovans.raw",
                               delim = " ")

geno_cage_bovans_load <- geno_cage_bovans[match(cage_bovans_load$data$X1, geno_cage_bovans$FID),]
assert_that(identical(geno_cage_bovans_load$FID, cage_bovans_load$data$X1))

geno_cage_bovans_weight <- geno_cage_bovans[match(cage_bovans_weight$data$X1, geno_cage_bovans$FID),]
assert_that(identical(geno_cage_bovans_weight$FID, cage_bovans_weight$data$X1))


geno_all_bovans <- read_delim("gwas/all_bovans.raw",
                              delim = " ")

geno_all_bovans_load <- geno_all_bovans[match(all_bovans_load$data$X1, geno_all_bovans$FID),]
assert_that(identical(geno_all_bovans_load$FID, all_bovans_load$data$X1))

geno_all_bovans_weight <- geno_all_bovans[match(all_bovans_weight$data$X1, geno_all_bovans$FID),]
assert_that(identical(geno_all_bovans_weight$FID, all_bovans_weight$data$X1))


geno_pen_lsl <- read_delim("gwas/pen_lsl.raw",
                              delim = " ")

geno_pen_lsl_load <- geno_pen_lsl[match(pen_lsl_load$data$X1, geno_pen_lsl$FID),]
assert_that(identical(geno_pen_lsl_load$FID, pen_lsl_load$data$X1))

geno_pen_lsl_weight <- geno_pen_lsl[match(pen_lsl_weight$data$X1, geno_pen_lsl$FID),]
assert_that(identical(geno_pen_lsl_weight$FID, pen_lsl_weight$data$X1))


geno_cage_lsl <- read_delim("gwas/cage_lsl.raw",
                               delim = " ")

geno_cage_lsl_load <- geno_cage_lsl[match(cage_lsl_load$data$X1, geno_cage_lsl$FID),]
assert_that(identical(geno_cage_lsl_load$FID, cage_lsl_load$data$X1))

geno_cage_lsl_weight <- geno_cage_lsl[match(cage_lsl_weight$data$X1, geno_cage_lsl$FID),]
assert_that(identical(geno_cage_lsl_weight$FID, cage_lsl_weight$data$X1))


geno_all_lsl <- read_delim("gwas/all_lsl.raw",
                           delim = " ")

geno_all_lsl_load <- geno_all_lsl[match(all_lsl_load$data$X1, geno_all_lsl$FID),]
assert_that(identical(geno_all_lsl_load$FID, all_lsl_load$data$X1))

geno_all_lsl_weight <- geno_all_lsl[match(all_lsl_weight$data$X1, geno_all_lsl$FID),]
assert_that(identical(geno_all_lsl_weight$FID, all_lsl_weight$data$X1))


## Pruned genotypes

snps_pruned_pen_bovans_load <- prune_snp_matrix(geno_pen_bovans_load[, -(1:6)])
snps_pruned_pen_bovans_weight <- prune_snp_matrix(geno_pen_bovans_weight[, -(1:6)])

snps_pruned_cage_bovans_load <- prune_snp_matrix(geno_cage_bovans_load[, -(1:6)])
snps_pruned_cage_bovans_weight <- prune_snp_matrix(geno_cage_bovans_weight[, -(1:6)])

snps_pruned_all_bovans_load <- prune_snp_matrix(geno_all_bovans_load[, -(1:6)])
snps_pruned_all_bovans_weight <- prune_snp_matrix(geno_all_bovans_weight[, -(1:6)])


snps_pruned_pen_lsl_load <- prune_snp_matrix(geno_pen_lsl_load[, -(1:6)])
snps_pruned_pen_lsl_weight <- prune_snp_matrix(geno_pen_lsl_weight[, -(1:6)])

snps_pruned_cage_lsl_load <- prune_snp_matrix(geno_cage_lsl_load[, -(1:6)])
snps_pruned_cage_lsl_weight <- prune_snp_matrix(geno_cage_lsl_weight[, -(1:6)])

snps_pruned_all_lsl_load <- prune_snp_matrix(geno_all_lsl_load[, -(1:6)])
snps_pruned_all_lsl_weight <- prune_snp_matrix(geno_all_lsl_weight[, -(1:6)])
