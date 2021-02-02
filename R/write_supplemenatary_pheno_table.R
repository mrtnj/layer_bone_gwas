library(dplyr)

pheno <- readRDS("outputs/pheno.Rds")

pheno <- filter(pheno, !is.na(cage.pen))

pheno$feed <- ifelse(grepl("kontroll", pheno$feed),
                     "control",
                     "treatment")

to_exclude <- c(grep("^ftir_", colnames(pheno), value = TRUE),
                "comb_g", "rna_sampled", "day")

pheno <- pheno[! colnames(pheno) %in% to_exclude]

write.csv(pheno,
          file = "tables/supplementary_data3_phenotypes.csv",
          quote = FALSE,
          row.names = FALSE)