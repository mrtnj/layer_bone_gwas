
library(assertthat)
library(dplyr)
library(readr)
library(tidyr)


## Prepare files for MTG2 and GTCA


## Fam file

all_load <- read_delim("gwas/all_load_adj/all_load_adj.fam",
                       delim = " ",
                       col_names = FALSE)

all_missing <- is.na(all_load$X6)

all_load <- all_load[!all_missing,]


write.table(all_load,
            file = "genomic_correlation/all.fam",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)


## GRM file

all_GRM <- as.data.frame(read_tsv("gwas/output/all_grm.sXX.txt", col_names = FALSE))

all_GRM <- all_GRM[!all_missing, !all_missing]

all_GRM$row_number <- 1:nrow(all_GRM)

long_GRM <- pivot_longer(all_GRM, -row_number, names_to = "col_number")

long_GRM$col_number <- sub(long_GRM$col_number, pattern = "X", replacement = "")


write.table(long_GRM,
            file = "genomic_correlation/grm.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)



## Phenotypes and covariates

pheno <- readRDS("outputs/pheno.Rds")

all_covar <- filter(pheno,
                    animal_id %in% all_load$X1 &
                    !is.na(load_N))[, c("animal_id", "group", "weight",
                                        "breed", "cage.pen", "load_N")]

all_covar <- all_covar[order(all_covar$animal_id),]

assert_that(identical(all_covar$animal_id, all_load$X1))


## Trait file

all_covar$pen_load <- ifelse(all_covar$cage.pen == "PEN",
                             all_load$X6,
                             NA)

all_covar$cage_load <- ifelse(all_covar$cage.pen == "CAGE",
                              all_load$X6,
                              NA)



write.table(all_covar[, c("animal_id", "animal_id",
                          "cage_load", "pen_load")],
            file = "genomic_correlation/pheno.txt",
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE)


## Categorical covariates

write.table(all_covar[, c("animal_id", "animal_id",
                          "breed")],
            file = "genomic_correlation/cc.txt",
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE)

## Quantitative covariates

write.table(all_covar[, c("animal_id", "animal_id",
                          "weight")],
            file = "genomic_correlation/qc.txt",
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE)
