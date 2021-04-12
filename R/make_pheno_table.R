library(dplyr)
library(readr)

## Collate phenotypes from all sources into one file


read_data <- function(filename) {
    read.table(filename,
               header = TRUE,
               sep = "\t")
}

## Breaking strength

breaking <- read_data("data/breaking_strength.txt")
colnames(breaking)[2] <- "load_N"


## Covariates

covariates <- lapply(c("data/cage_or_pen.txt",
                       "data/feeds.txt",
                       "data/breeds.txt",
                       "data/groups.txt"),
                     read_data)

covariates <- Reduce(full_join, covariates)


## Weight

weight <- read_data("data/weights.txt")



## pQCT data

ct_distal <- read_data("data/micro_distal_avg.txt")

colnames(ct_distal)[5:58] <- paste("ct_distal_",
                                   colnames(ct_distal)[5:58],
                                   sep = "")

ct_mid <- read_data("data/micro_midshaft_avg.txt")

colnames(ct_mid)[5:58] <- paste("ct_mid_",
                                colnames(ct_mid)[5:58],
                                sep = "")

ct <- inner_join(ct_distal[, -(2:4)],
                 ct_mid[, -(2:4)])


## pQCT principal components

ct_trait_names <- c("TOT_CNT", "TOT_DEN",
                    "TRAB_CNT", "TRAB_DEN",
                    "CRT_CNT", "CRT_DEN",
                    "CRT_THK_C", "OBJECTLEN")

ct_ix <- which(colnames(ct) %in%
                   c(paste("ct_mid_", ct_trait_names, sep = ""),
                     paste("ct_distal_", ct_trait_names, sep = "")))


ct_animals <- na.exclude(ct[, c(1, ct_ix)])$animal_id

ct_complete <- filter(ct, animal_id %in% ct_animals)

ct_pca <- prcomp(ct_complete[, ct_ix], scale = TRUE)

ct_pca_data <- data.frame(animal_id = ct_complete$animal_id,
                          ct_pc1 = ct_pca$x[, 1],
                          ct_pc2 = ct_pca$x[, 2],
                          ct_pc3 = ct_pca$x[, 3])



## TGA data

tga <- read_tsv("data/tga.txt",
                na = c("", "--"))
colnames(tga) <- sub(colnames(tga), pattern = "/", replacement = "_over_")


pheno <- Reduce(full_join, list(breaking, covariates, weight,
                                ct,
                                ct_pca_data,
                                tga))


## Order phenotype data
pheno <- pheno[order(pheno$animal_id),]

saveRDS(pheno,
        file = "outputs/pheno.Rds")



pheno <- filter(pheno, !is.na(cage.pen))

write.table(pheno,
            file = "tables/supplementary_data3_phenotypes.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)