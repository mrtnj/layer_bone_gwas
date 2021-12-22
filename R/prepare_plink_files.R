
## Read Access phenotype and SNP genotype files to create PLINK files

library(assertthat)
library(dplyr)
library(readr)

source("R/preparation_helper_functions.R")

pheno <- readRDS("outputs/pheno.Rds")


## Fam and covariate files for cage/pen-separated GWAS

pen <- filter(pheno,
              cage.pen == "PEN" &
              !is.na(breed))

cage <- filter(pheno,
               cage.pen == "CAGE" &
               !is.na(breed))

pen <- filter(pheno,
              cage.pen == "PEN" &
                  !is.na(breed))

bovans <- filter(pheno,
                 !is.na(cage.pen) &
                     breed == "Bovans")

LSL <- filter(pheno,
              !is.na(cage.pen) &
                  breed == "LSL")

pen_bovans <- filter(pheno,
                     cage.pen == "PEN" &
                         breed == "Bovans")

pen_LSL <- filter(pheno,
                  cage.pen == "PEN" &
                      breed == "LSL")

cage_bovans <- filter(pheno,
                      cage.pen == "CAGE" &
                          breed == "Bovans")

cage_LSL <- filter(pheno,
                   cage.pen == "CAGE" &
                       breed == "LSL")




geno <- read_tsv("data/RA-1698_181010_SNPGenotypeExport_PCF_TOP.txt",
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


saveRDS(geno_pruned,
        file = "outputs/geno.Rds")


## Format genotypes for converting to plink compound ped

geno_compound <- as.data.frame(geno_pruned)

for (col_ix in 2:ncol(geno)) {
    geno_compound[, col_ix] <- sub(geno_compound[, col_ix],
                                   pattern = "/",
                                   replacement = " ")
}

## Genotype of all

geno_all <- geno_compound[na.exclude(match(pheno$animal_id, geno_compound$individual)),]

saveRDS(geno_all,
        file = "outputs/geno_all.Rds")


geno_all_bovans <- geno_compound[na.exclude(match(bovans$animal_id, geno_compound$individual)),]
geno_all_LSL <- geno_compound[na.exclude(match(LSL$animal_id, geno_compound$individual)),]



## Subset genotypes based on pen/cage

geno_pen <- geno_compound[na.exclude(match(pen$animal_id, geno_compound$individual)),]
geno_cage <- geno_compound[na.exclude(match(cage$animal_id, geno_compound$individual)),]

geno_pen_bovans <- geno_compound[na.exclude(match(pen_bovans$animal_id, geno_compound$individual)),]
geno_cage_bovans <- geno_compound[na.exclude(match(cage_bovans$animal_id, geno_compound$individual)),]

geno_pen_LSL <- geno_compound[na.exclude(match(pen_LSL$animal_id, geno_compound$individual)),]
geno_cage_LSL <- geno_compound[na.exclude(match(cage_LSL$animal_id, geno_compound$individual)),]



## Subset phenotypes based on genotypes

all_genotyped <- filter(pheno,
                        animal_id %in% geno_all$individual)

pen_genotyped <- filter(pen,
                        animal_id %in% geno_pen$individual)
cage_genotyped <- filter(cage,
                         animal_id %in% geno_cage$individual)


all_bovans_genotyped <- filter(bovans,
                               animal_id %in% geno_all_bovans$individual)

all_LSL_genotyped <- filter(LSL,
                            animal_id %in% geno_all_LSL$individual)


pen_bovans_genotyped <- filter(pen_bovans,
                               animal_id %in% geno_pen_bovans$individual)
cage_bovans_genotyped <- filter(cage_bovans,
                                animal_id %in% geno_cage_bovans$individual)


pen_LSL_genotyped <- filter(pen_LSL,
                            animal_id %in% geno_pen_LSL$individual)
cage_LSL_genotyped <- filter(cage_LSL,
                                animal_id %in% geno_cage_LSL$individual)



## Check ids
assert_that(identical(all_genotyped$animal_id,
                      geno_all$individual))

assert_that(identical(pen_genotyped$animal_id,
                      geno_pen$individual))
assert_that(identical(cage_genotyped$animal_id,
                      geno_cage$individual))

assert_that(identical(all_bovans_genotyped$animal_id,
                      geno_all_bovans$individual))
assert_that(identical(all_LSL_genotyped$animal_id,
                      geno_all_LSL$individual))


assert_that(identical(pen_bovans_genotyped$animal_id,
                      geno_pen_bovans$individual))
assert_that(identical(cage_bovans_genotyped$animal_id,
                      geno_cage_bovans$individual))

assert_that(identical(pen_LSL_genotyped$animal_id,
                      geno_pen_LSL$individual))
assert_that(identical(cage_LSL_genotyped$animal_id,
                      geno_cage_LSL$individual))


## Traits of interest

traits <- c("load_N", "weight", 
            "ct_pc1", "ct_pc2", "ct_pc3",
            "WaterLost_CB", "OMLost_CB", "CO2Lost_CB",
            "Phosphates_CB", "Mineral_CB",
            "Phosphates_over_OM_CB", "CO3_over_Phosphates_CB",
            "WaterLost_MB", "OMLost_MB",
            "CO2Lost_MB", "Phosphates_MB",
            "Mineral_MB", "Phosphates_over_OM_MB",
            "CO3_over_Phosphates_MB")


## Create fam files


for (trait_ix in 1:length(traits)) {
 
    fam_all <- fam(all_genotyped, traits[trait_ix])
    write_plink(fam_all,
                paste("gwas/fam_all_", traits[trait_ix], ".fam", sep = ""))
    
    fam_pen <- fam(pen_genotyped, traits[trait_ix])
    write_plink(fam_pen,
                paste("gwas/fam_pen_", traits[trait_ix], ".fam", sep = ""))
    
    fam_cage <- fam(cage_genotyped, traits[trait_ix])
    write_plink(fam_cage,
                paste("gwas/fam_cage_", traits[trait_ix], ".fam", sep = ""))
    
    assert_that(identical(geno_all$individual, fam_all$fid))
    assert_that(identical(geno_pen$individual, fam_pen$fid))
    assert_that(identical(geno_cage$individual, fam_cage$fid))
    
}

traits_breed <- traits[1:2]

for (trait_ix in 1:length(traits_breed)) {
    
    fam_pen_bovans <- fam(pen_bovans_genotyped, traits[trait_ix])
    write_plink(fam_pen_bovans,
                paste("gwas/fam_pen_bovans_", traits[trait_ix], ".fam", sep = ""))
    
    fam_cage_bovans <- fam(cage_bovans_genotyped, traits[trait_ix])
    write_plink(fam_cage_bovans,
                paste("gwas/fam_cage_bovans_", traits[trait_ix], ".fam", sep = ""))
    
    fam_pen_LSL <- fam(pen_LSL_genotyped, traits[trait_ix])
    write_plink(fam_pen_bovans,
                paste("gwas/fam_pen_lsl_", traits[trait_ix], ".fam", sep = ""))
    
    fam_cage_LSL <- fam(cage_LSL_genotyped, traits[trait_ix])
    write_plink(fam_cage_LSL,
                paste("gwas/fam_cage_lsl_", traits[trait_ix], ".fam", sep = ""))
    
}



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

ped_all_bovans <- ped(geno_all_bovans)
ped_all_LSL <- ped(geno_all_LSL)

ped_pen <- ped(geno_pen)
ped_cage <- ped(geno_cage)

ped_pen_bovans <- ped(geno_pen_bovans)
ped_cage_bovans <- ped(geno_cage_bovans)

ped_pen_LSL <- ped(geno_pen_LSL)
ped_cage_LSL <- ped(geno_cage_LSL)


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



## Write out everything

write_ped <- function(x, filename) write.table(x,
                                               file = filename,
                                               quote = FALSE,
                                               row.names = FALSE,
                                               col.names = FALSE,
                                               na = "0 0")

## Covariates

write_plink(covar_all_cagepen_breed_weight, "gwas/covar_all_breed_weight.txt")
write_plink(covar_all_cagepen_breed, "gwas/covar_all_breed.txt")

write_plink(covar_pen_breed_weight, "gwas/covar_pen_breed_weight.txt")
write_plink(covar_cage_breed_weight, "gwas/covar_cage_breed_weight.txt")
write_plink(covar_pen_breed, "gwas/covar_pen_breed.txt")
write_plink(covar_cage_breed, "gwas/covar_cage_breed.txt")
write_plink(covar_pen_weight, "gwas/covar_pen_weight.txt")
write_plink(covar_cage_weight, "gwas/covar_cage_weight.txt")


## Ped

write_ped(ped_all, "gwas/all.ped")

write_ped(ped_all_bovans, "gwas/all_bovans.ped")
write_ped(ped_all_LSL, "gwas/all_lslL.ped")

write_ped(ped_pen, "gwas/pen.ped")
write_ped(ped_cage, "gwas/cage.ped")

write_ped(ped_pen_bovans, "gwas/pen_bovans.ped")
write_ped(ped_cage_bovans, "gwas/cage_bovans.ped")

write_ped(ped_pen_LSL, "gwas/pen_lsl.ped")
write_ped(ped_cage_LSL, "gwas/cage_lsl.ped")

