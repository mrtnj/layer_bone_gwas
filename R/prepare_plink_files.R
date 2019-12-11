
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



## Fam and covariate files for cage/pen-separated and breed/cage/pen-separated GWAS

pen <- filter(pheno,
              cage.pen == "PEN" &
              !is.na(breed))

cage <- filter(pheno,
               cage.pen == "CAGE" &
               !is.na(breed))

bovans_pen <- filter(pheno,
                     cage.pen == "PEN" &
                     breed == "Bovans")

bovans_cage <- filter(pheno,
                     cage.pen == "CAGE" &
                     breed == "Bovans")

lsl_pen <- filter(pheno,
                  cage.pen == "PEN" &
                  breed == "LSL")

lsl_cage <- filter(pheno,
                   cage.pen == "CAGE" &
                   breed == "LSL")


geno <- read_tsv("RA-1698_181010_ResultReport/RA-1698_181010_ResultReport_PCF_TOP/RA-1698_181010_SNPGenotypeExport_PCF_TOP.txt",
                 col_types = cols(.default = "c", individual = "n"))

geno <- geno[order(geno$individual),]



## Get the map

map <- read_delim("gwas/map.map",
                  delim = " ",
                  col_names = FALSE,
                  col_types = "cccn")

col_ix <- c(1, na.exclude(match(map$X2, colnames(geno))))

geno <- geno[, col_ix]

assert_that(identical(colnames(geno)[-1],
            map$X2))


## High missingness individuals to exclude

ids_high_missingness <- scan("outputs/ids_high_missingness.txt")

ids_strange_breed <- scan("outputs/ids_strange_breed_assignment.txt")

geno_pruned <- filter(geno, !(individual %in% c(ids_high_missingness,
                                                ids_strange_breed)))


## Format genotypes for converting to plink compound ped

geno_compound <- as.data.frame(geno_pruned)

for (col_ix in 2:ncol(geno)) {
    geno_compound[, col_ix] <- sub(geno_compound[, col_ix],
                                   pattern = "/",
                                   replacement = " ")
}

## Genotype of all

geno_all <- geno_compound[na.exclude(match(pheno$animal_id, geno_compound$individual)),]


## Subset genotypes based on pen/cage

geno_pen <- geno_compound[na.exclude(match(pen$animal_id, geno_compound$individual)),]
geno_cage <- geno_compound[na.exclude(match(cage$animal_id, geno_compound$individual)),]

geno_bovans_pen <- geno_compound[na.exclude(match(bovans_pen$animal_id,
                                                  geno_compound$individual)),]
geno_bovans_cage <- geno_compound[na.exclude(match(bovans_cage$animal_id,
                                                   geno_compound$individual)),]

geno_lsl_pen <- geno_compound[na.exclude(match(lsl_pen$animal_id,
                                               geno_compound$individual)),]
geno_lsl_cage <- geno_compound[na.exclude(match(lsl_cage$animal_id,
                                                geno_compound$individual)),]


## Subset phenotypes based on genotypes

all_genotyped <- filter(pheno,
                        animal_id %in% geno_all$individual)

pen_genotyped <- filter(pen,
                        animal_id %in% geno_pen$individual)
cage_genotyped <- filter(cage,
                         animal_id %in% geno_cage$individual)

bovans_pen_genotyped <- filter(bovans_pen,
                               animal_id %in% geno_bovans_pen$individual)
bovans_cage_genotyped <- filter(bovans_cage,
                                animal_id %in% geno_bovans_cage$individual)

lsl_pen_genotyped <- filter(lsl_pen,
                            animal_id %in% geno_lsl_pen$individual)
lsl_cage_genotyped <- filter(lsl_cage,
                             animal_id %in% geno_lsl_cage$individual)



## Check ids
assert_that(identical(all_genotyped$animal_id,
                      geno_all$individual))

assert_that(identical(pen_genotyped$animal_id,
                      geno_pen$individual))
assert_that(identical(cage_genotyped$animal_id,
                      geno_cage$individual))

assert_that(identical(bovans_pen_genotyped$animal_id,
                      geno_bovans_pen$individual))
assert_that(identical(bovans_cage_genotyped$animal_id,
                      geno_bovans_cage$individual))

assert_that(identical(lsl_pen_genotyped$animal_id,
                      geno_lsl_pen$individual))
assert_that(identical(lsl_cage_genotyped$animal_id,
                      geno_lsl_cage$individual))


## Create fam files

fam <- function(df, trait) data.frame(fid = df$animal_id,
                                      iid = df$animal_id,
                                      father = 0,
                                      mother = 0,
                                      sex = 2,
                                      as.data.frame(df[, trait]))

fam_all_load <- fam(all_genotyped, "load_N")
fam_all_weight <- fam(all_genotyped, "weight")
fam_all_comb <- fam(all_genotyped, "comb_g")

fam_pen_load <- fam(pen_genotyped, "load_N")
fam_pen_weight <- fam(pen_genotyped, "weight")
fam_pen_comb <- fam(pen_genotyped, "comb_g")

fam_cage_load <- fam(cage_genotyped, "load_N")
fam_cage_weight <- fam(cage_genotyped, "weight")
fam_cage_comb <- fam(cage_genotyped, "comb_g")

fam_bovans_pen_load <- fam(bovans_pen_genotyped, "load_N")
fam_bovans_pen_weight <- fam(bovans_pen_genotyped, "weight")
fam_bovans_pen_comb <- fam(bovans_pen_genotyped, "comb_g")

fam_lsl_pen_load <- fam(lsl_pen_genotyped, "load_N")
fam_lsl_pen_weight <- fam(lsl_pen_genotyped, "weight")
fam_lsl_pen_comb <- fam(lsl_pen_genotyped, "comb_g")

fam_bovans_cage_load <- fam(bovans_cage_genotyped, "load_N")
fam_bovans_cage_weight <- fam(bovans_cage_genotyped, "weight")
fam_bovans_cage_comb <- fam(bovans_cage_genotyped, "comb_g")

fam_lsl_cage_load <- fam(lsl_cage_genotyped, "load_N")
fam_lsl_cage_weight <- fam(lsl_cage_genotyped, "weight")
fam_lsl_cage_comb <- fam(lsl_cage_genotyped, "comb_g")



## Create ped files

ped <- function(df) data.frame(fid = df$individual,
                               iid = df$individual,
                               father = 0,
                               mother = 0,
                               sex = 2,
                               trait = 1,
                               df[, -1],
                               stringsAsFactors = FALSE)

ped_all <- ped(geno_all)

ped_pen <- ped(geno_pen)
ped_cage <- ped(geno_cage)

ped_bovans_pen <- ped(geno_bovans_pen)
ped_bovans_cage <- ped(geno_bovans_cage)

ped_lsl_pen <- ped(geno_lsl_pen)
ped_lsl_cage <- ped(geno_lsl_cage)


## Create covariate tables

covar_all_cagepen_breed_weight <- model.matrix(~ cage.pen + breed + weight,
                                               data = all_genotyped)

covar_all_cagepen_breed <- model.matrix(~ cage.pen + breed,
                                        data = all_genotyped)


covar_pen_breed_weight <- model.matrix(~ breed + weight,
                                       data = pen_genotyped)

covar_pen_breed <- model.matrix(~ breed,
                                data = pen_genotyped)

covar_pen_weight <- model.matrix(~ weight,
                                data = pen_genotyped)

covar_cage_breed_weight <- model.matrix(~ breed + weight,
                                        data = cage_genotyped)

covar_cage_breed <- model.matrix(~ breed,
                                 data = cage_genotyped)

covar_cage_weight <- model.matrix(~ weight,
                                 data = cage_genotyped)

covar_bovans_pen_weight <- model.matrix(~ weight,
                                        data = bovans_pen_genotyped)

covar_bovans_cage_weight <- model.matrix(~ weight,
                                         data = bovans_cage_genotyped)

covar_lsl_pen_weight <- model.matrix(~ weight,
                                     data = lsl_pen_genotyped)

covar_lsl_cage_weight <- model.matrix(~ weight,
                                      data = lsl_cage_genotyped)


## Write out everything

write_plink <- function(x, filename) write.table(x,
                                                 file = filename,
                                                 quote = FALSE,
                                                 row.names = FALSE,
                                                 col.names = FALSE)

write_ped <- function(x, filename) write.table(x,
                                               file = filename,
                                               quote = FALSE,
                                               row.names = FALSE,
                                               col.names = FALSE,
                                               na = "0 0")
## Fam

write_plink(fam_all_load, "gwas/fam_all_load.fam")
write_plink(fam_all_comb, "gwas/fam_all_comb.fam")
write_plink(fam_all_weight, "gwas/fam_all_weight.fam")

write_plink(fam_pen_load, "gwas/fam_pen_load.fam")
write_plink(fam_cage_load, "gwas/fam_cage_load.fam")
write_plink(fam_pen_weight, "gwas/fam_pen_weight.fam")
write_plink(fam_cage_weight, "gwas/fam_cage_weight.fam")
write_plink(fam_pen_comb, "gwas/fam_pen_comb.fam")
write_plink(fam_cage_comb, "gwas/fam_cage_comb.fam")

write_plink(fam_bovans_pen_load, "gwas/fam_bovans_pen_load.fam")
write_plink(fam_bovans_pen_weight, "gwas/fam_bovans_pen_weight.fam")
write_plink(fam_bovans_pen_comb, "gwas/fam_bovans_pen_comb.fam")

write_plink(fam_lsl_pen_load, "gwas/fam_lsl_pen_load.fam")
write_plink(fam_lsl_pen_weight, "gwas/fam_lsl_pen_weight.fam")
write_plink(fam_lsl_pen_comb, "gwas/fam_lsl_pen_comb.fam")

write_plink(fam_bovans_cage_load, "gwas/fam_bovans_cage_load.fam")
write_plink(fam_bovans_cage_weight, "gwas/fam_bovans_cage_weight.fam")
write_plink(fam_bovans_cage_comb, "gwas/fam_bovans_cage_comb.fam")

write_plink(fam_lsl_cage_load, "gwas/fam_lsl_cage_load.fam")
write_plink(fam_lsl_cage_weight, "gwas/fam_lsl_cage_weight.fam")
write_plink(fam_lsl_cage_comb, "gwas/fam_lsl_cage_comb.fam")


## Covariates

write_plink(covar_all_cagepen_breed_weight, "gwas/covar_all_cagepen_breed_weight.txt")
write_plink(covar_all_cagepen_breed, "gwas/covar_all_cagepen_breed.txt")

write_plink(covar_pen_breed_weight, "gwas/covar_pen_breed_weight.txt")
write_plink(covar_cage_breed_weight, "gwas/covar_cage_breed_weight.txt")
write_plink(covar_pen_breed, "gwas/covar_pen_breed.txt")
write_plink(covar_cage_breed, "gwas/covar_cage_breed.txt")
write_plink(covar_pen_weight, "gwas/covar_pen_weight.txt")
write_plink(covar_cage_weight, "gwas/covar_cage_weight.txt")

write_plink(covar_bovans_pen_weight, "gwas/covar_bovans_pen_weight.txt")
write_plink(covar_lsl_pen_weight, "gwas/covar_lsl_pen_weight.txt")
write_plink(covar_bovans_cage_weight, "gwas/covar_bovans_cage_weight.txt")
write_plink(covar_lsl_cage_weight, "gwas/covar_lsl_cage_weight.txt")


## Ped

write_ped(ped_all, "gwas/all.ped")

write_ped(ped_pen, "gwas/pen.ped")
write_ped(ped_cage, "gwas/cage.ped")

write_ped(ped_bovans_pen, "gwas/bovans_pen.ped")
write_ped(ped_bovans_cage, "gwas/bovans_cage.ped")
write_ped(ped_lsl_pen, "gwas/lsl_pen.ped")
write_ped(ped_lsl_cage, "gwas/lsl_cage.ped")
