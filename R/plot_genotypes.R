

## Exploratory plots of genotypes

library(dplyr)
library(ggplot2)
library(patchwork)
library(purrr)
library(readr)
library(tidyr)

geno <- read_tsv("data/RA-1698_181010_SNPGenotypeExport_PCF_TOP.txt",
                 col_types = cols(.default = "c", individual = "n"))


## Recode genotypes for PCA

recode_locus <- function(loc) {
    if (is.na(loc[1])) {
        r <- NA
    } else if (loc[1] == loc[2] & loc[1] == alleles[1]) {
        r <- 0
    } else if (loc[1] == loc[2] & loc[1] == alleles[2]) {
        r <- 2
    } else if (loc[1] != loc[2]) {
        r <-  1
    }
    r
}

recoded <- as.data.frame(geno)
monomorphic <- logical(ncol(geno))

for (col_ix in 2:ncol(geno)) {
    geno_split <- strsplit(recoded[, col_ix],
                           split = "/")
    alleles <- na.exclude(unique(unlist(geno_split)))
    
    if (length(alleles) == 1) {
        ## Monomorphic
        monomorphic[col_ix] <- TRUE
    } else {
        recoded[, col_ix] <- unlist(lapply(geno_split, recode_locus))
    }
}

## Missingness

missing <- numeric(nrow(recoded))

for (row_ix in 1:nrow(recoded)) {
    missing[row_ix] <- sum(is.na(recoded[row_ix,]))
}

missing_df <- data.frame(animal_id = recoded$individual,
                         missing = missing)


## Restrict to segregating sites

segregating <- recoded[, !monomorphic]


## Frequencies

f <- colSums(segregating[,-1], na.rm = TRUE)/2/nrow(segregating)
maf <- ifelse(f > 0.5, 1-f, f)

plot_maf <- qplot(maf)


## Stupid imputation by drawing from allele frequency

sample_from_f <- function(n, f) {
    rbinom(n = n,
           size = 2,
           prob = f)
}

imputed <- segregating

for (col_ix in 2:ncol(imputed)) {

    ix_to_fill <- which(is.na(imputed[, col_ix]))

    if (length(ix_to_fill) > 0) {
        imputed[ix_to_fill, col_ix] <-
            sample_from_f(length(ix_to_fill),
                          f[col_ix - 1])
    }
}


## PCA and plots

pca <- prcomp(imputed[,-1])

pc_data <- data.frame(animal_id = geno$individual, pca$x)

pheno <- readRDS("outputs/pheno.Rds")

pheno <- filter(pheno, !is.na(cage.pen))


pc_pheno <- Reduce(function(x, y) inner_join(x, y, by = "animal_id"),
                   list(pheno, pc_data, missing_df))


plot_pc12 <- qplot(x = PC1, y = PC2, colour = breed, data = pc_pheno)

plot_zoomed <- plot_pc12 + coord_cartesian(ylim = c(-5, 5))

plot_pc_missing <- qplot(x = PC1, y = PC2, colour = missing, data = pc_pheno)

plot_zoomed_missing <- plot_pc_missing + coord_cartesian(ylim = c(-5, 5))

plot_breed_missing <- qplot(x = PC1, y = missing, colour = breed, data = pc_pheno) + ylim(0, 1000)


## Find animals that look extreme or misplaced based on their PCs

ids_extreme <- pc_pheno$animal_id[pc_pheno$PC2 < -10]
ids_strange_bovans <- pc_pheno$animal_id[pc_pheno$breed == "Bovans" &
                                         pc_pheno$PC1 > 0]
ids_strange_lsl <- pc_pheno$animal_id[pc_pheno$breed == "LSL" &
                                      pc_pheno$PC1 < 0]

ids_complete <- pc_pheno$animal_id[pc_pheno$missing < 1000]

ids_missingness <- pc_pheno$animal_id[pc_pheno$missing > 1000]

write(ids_extreme,
      file = "outputs/ids_extreme_pcs.txt",
      sep = "\n")

write(c(ids_strange_bovans, ids_strange_lsl),
      file = "outputs/ids_strange_breed_assignment.txt",
      sep = "\n")

write(ids_missingness,
      file = "outputs/ids_high_missingness.txt",
      sep = "\n")


## Clean up missingness

imputed_pruned <- filter(imputed, individual %in% ids_complete)


pca_pruned <- prcomp(imputed_pruned[,-1])

pc_data_pruned <- data.frame(animal_id = imputed_pruned$individual, pca_pruned$x)

pc_pheno_pruned <- Reduce(function(x, y) inner_join(x, y, by = "animal_id"),
                          list(pheno, pc_data_pruned, missing_df))

plot_pc12_pruned <- qplot(x = PC1, y = PC2, colour = breed, data = pc_pheno_pruned) +
    scale_colour_manual(values = c("blue", "red"),
                        name = "") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    ggtitle("Principal components of genotypes")
                


long_pc <- pivot_longer(pc_pheno_pruned,
                        paste("PC", 1:50, sep = ""),
                        names_to = "component")

plot_pcs_breed <- qplot(x = component, y = value, colour = breed, data = long_pc)
plot_pcs_pen <- qplot(x = component, y = value, colour = pen.cage, data = long_pc)




pdf("figures/plot_pcs.pdf")
print(plot_pc12_pruned)
dev.off()



## Between breed comparison of allele frequency

LSL_ids <- na.exclude(pheno$animal_id[pheno$breed == "LSL"])
bovans_ids <- na.exclude(pheno$animal_id[pheno$breed == "Bovans"])


geno_LSL <- filter(segregating, individual %in% LSL_ids)

geno_bovans <- filter(segregating, individual %in% bovans_ids)



nonmissing_LSL <- map_dbl(geno_LSL[, -1], function(g) sum(!is.na(g)))

nonmissing_bovans <- map_dbl(geno_bovans[, -1], function(g) sum(!is.na(g)))


f_LSL <- colSums(geno_LSL[, -1], na.rm = TRUE)/2/nonmissing_LSL

f_bovans <- colSums(geno_bovans[, -1], na.rm = TRUE)/2/nonmissing_bovans


plot_frequency_comparion <- qplot(x = f_bovans,
                                  y = f_LSL) /
  qplot(x = f_bovans - f_LSL)



cor(f_bovans, f_LSL, use = "p")
