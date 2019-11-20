
## Read Access dump files and create phenotype files

library(dplyr)


read_data <- function(filename) {
    read.table(filename,
               header = TRUE,
               dec = ",",
               sep = "\t")
}

## Breaking strenght

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
