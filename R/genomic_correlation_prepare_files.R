
library(assertthat)
library(dplyr)
library(readr)
library(tidyr)


## Prepare files for GTCA

system("mkdir genomic_correlation")
system("mkdir genomic_correlation/load_N")
system("mkdir genomic_correlation/weight")
system("mkdir genomic_correlation/comb")


## Fam file

all_load <- read_delim("gwas/all_load_N/all_load_N.fam",
                       delim = " ",
                       col_names = FALSE)

all_missing <- is.na(all_load$X6)

##all_load <- all_load[!all_missing,]



write.table(all_load,
            file = "genomic_correlation/load_N/all.fam",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)


all_weight <- read_delim("gwas/all_weight/all_weight.fam",
                       delim = " ",
                       col_names = FALSE)

all_weight_missing <- is.na(all_weight$X6)

## all_weight <- all_weight[!all_weight_missing,]


write.table(all_weight,
            file = "genomic_correlation/weight/all.fam",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)





## Phenotypes and covariates

pheno <- readRDS("outputs/pheno.Rds")

load_covar <- filter(pheno,
                     animal_id %in% all_load$X1)[, c("animal_id", "group", "weight",
                                                     "breed", "cage.pen", "load_N", "comb_g")]

assert_that(identical(load_covar$animal_id, all_load$X1))


weight_covar <- filter(pheno,
                       animal_id %in% all_weight$X1 &
                           !is.na(weight))[, c("animal_id", "group", "weight",
                                               "breed", "cage.pen")]

assert_that(identical(weight_covar$animal_id, all_weight$X1))



## Trait file

load_covar$pen_load <- ifelse(load_covar$cage.pen == "PEN",
                              load_covar$load_N,
                              NA)

load_covar$cage_load <- ifelse(load_covar$cage.pen == "CAGE",
                               load_covar$load_N,
                               NA)


load_covar$pen_comb <- ifelse(load_covar$cage.pen == "PEN",
                              load_covar$comb_g,
                              NA)

load_covar$cage_comb <- ifelse(load_covar$cage.pen == "CAGE",
                               load_covar$comb_g,
                               NA)


weight_covar$pen_weight <- ifelse(weight_covar$cage.pen == "PEN",
                                  weight_covar$weight,
                                  NA)

weight_covar$cage_weight <- ifelse(weight_covar$cage.pen == "CAGE",
                                   weight_covar$weight,
                                   NA)




write.table(load_covar[, c("animal_id", "animal_id",
                          "cage_load", "pen_load")],
            file = "genomic_correlation/load_N/pheno.txt",
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE)

write.table(load_covar[, c("animal_id", "animal_id",
                          "cage_comb", "pen_comb")],
            file = "genomic_correlation/comb/pheno_comb.txt",
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE)

write.table(weight_covar[, c("animal_id", "animal_id",
                             "cage_weight", "pen_weight")],
            file = "genomic_correlation/weight/pheno.txt",
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE)


## Categorical covariates

write.table(load_covar[, c("animal_id", "animal_id",
                          "breed")],
            file = "genomic_correlation/load_N/cc.txt",
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE)

## Quantitative covariates

write.table(load_covar[, c("animal_id", "animal_id",
                          "weight")],
            file = "genomic_correlation/load_N/qc.txt",
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE)


## Categorical covariates for weight

write.table(weight_covar[, c("animal_id", "animal_id",
                           "breed")],
            file = "genomic_correlation/weight/cc.txt",
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE)


