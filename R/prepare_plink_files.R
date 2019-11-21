
## Read Access dump files and SNP genotype files to create PLINK files for GEMMA

library(assertthat)
library(dplyr)
library(readr)


read_data <- function(filename) {
    read.table(filename,
               header = TRUE,
               dec = ",",
               sep = "\t")
}

## Breaking strength

breaking <- read_data("mdb_dump/breaking_strength.txt")
colnames(breaking)[2] <- "load_N"


## Covariates

covariates <- lapply(c("mdb_dump/cage_or_pen.txt",
                       "mdb_dump/feeds.txt",
                       "mdb_dump/breeds.txt",
                       "mdb_dump/groups.txt"),
                     read_data)

covariates <- Reduce(full_join, covariates)


## Weight

weight <- read_data("mdb_dump/weights.txt")


## Comb size

comb <- read_data("mdb_dump/combs.txt")


pheno <- Reduce(full_join, list(breaking, covariates, weight, comb))

saveRDS(pheno,
        file = "outputs/pheno.Rds")


## Order phenotype data
pheno <- pheno[order(pheno$animal_id),]


## Fam and covariate files for cage/pen-separated GWAS

pen <- filter(pheno,
              cage.pen == "PEN" &
              !is.na(breed))

cage <- filter(pheno,
               cage.pen == "CAGE" &
               !is.na(breed))



geno <- read_tsv("RA-1698_181010_ResultReport/RA-1698_181010_ResultReport_PCF_TOP/RA-1698_181010_SNPGenotypeExport_PCF_TOP.txt",
                 col_types = cols(.default = "c", individual = "n"))

geno <- geno[order(geno$individual),]


## High missingness individuals to exclude

ids_high_missingness <- scan("outputs/ids_high_missingness.txt")

geno_pruned <- filter(geno, !(individual %in% ids_high_missingness))


## Format genotypes for converting to plink compound ped

geno_compound <- as.data.frame(geno_pruned)

for (col_ix in 2:ncol(geno)) {
    geno_compound[, col_ix] <- sub(geno_compound[, col_ix],
                                   pattern = "/",
                                   replacement = "")
}


## Subset genotypes based on pen/cage

geno_pen <- geno_compound[na.exclude(match(pen$animal_id, geno_compound$individual)),]
geno_cage <- geno_compound[na.exclude(match(cage$animal_id, geno_compound$individual)),]



## Subset phenotypes based on genotypes

pen_genotyped <- filter(pen,
                        animal_id %in% geno_pen$individual)
cage_genotyped <- filter(cage,
                         animal_id %in% geno_cage$individual)


## Check ids
assert_that(identical(pen_genotyped$animal_id,
                      geno_pen$individual))
assert_that(identical(cage_genotyped$animal_id,
                      geno_cage$individual))


## Create fam files

fam <- function(df, trait) data.frame(fid = df$animal_id,
                                      iid = df$animal_id,
                                      father = 0,
                                      mother = 0,
                                      sex = 2,
                                      as.data.frame(df[, trait]))

fam_pen_load <- fam(pen_genotyped, "load_N")
fam_pen_weight <- fam(pen_genotyped, "weight")
fam_pen_comb <- fam(pen_genotyped, "comb_g")

fam_cage_load <- fam(cage_genotyped, "load_N")
fam_cage_weight <- fam(cage_genotyped, "weight")
fam_cage_comb <- fam(cage_genotyped, "comb_g")



## Create ped files

ped <- function(df) data.frame(fid = df$individual,
                               iid = df$individual,
                               father = 0,
                               mother = 0,
                               sex = 2,
                               df[, -1],
                               stringsAsFactors = FALSE)

ped_pen <- ped(geno_pen)
ped_cage <- ped(geno_cage)


## Create covariate tables

covar_pen_breed_weight <- model.matrix(~ breed + weight,
                                       data = pen_genotyped)

covar_pen_breed <- model.matrix(~ breed,
                                data = pen_genotyped)

covar_cage_breed_weight <- model.matrix(~ breed + weight,
                                        data = cage_genotyped)

covar_cage_breed <- model.matrix(~ breed,
                                data = cage_genotyped)


## Write out everything

write_plink <- function(x, filename) write.table(x,
                                                 file = filename,
                                                 quote = FALSE,
                                                 row.names = FALSE,
                                                 col.names = FALSE,
                                                 na = "-9")

write_plink(fam_pen_load, "gwas/fam_pen_load.fam")
write_plink(fam_cage_load, "gwas/fam_cage_load.fam")
write_plink(fam_pen_weight, "gwas/fam_pen_weight.fam")
write_plink(fam_cage_weight, "gwas/fam_cage_weight.fam")
write_plink(fam_pen_comb, "gwas/fam_pen_comb.fam")
write_plink(fam_cage_comb, "gwas/fam_cage_comb.fam")

write_plink(covar_pen_breed_weight, "gwas/covar_pen_breed_weight.txt")
write_plink(covar_cage_breed_weight, "gwas/covar_cage_breed_weight.txt")
write_plink(covar_pen_breed, "gwas/covar_pen_breed.txt")
write_plink(covar_cage_breed, "gwas/covar_cage_breed.txt")

write_plink(ped_pen, "gwas/pen.ped")
write_plink(ped_cage, "gwas/cage.ped")
