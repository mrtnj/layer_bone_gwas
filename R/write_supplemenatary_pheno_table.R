library(dplyr)

pheno <- readRDS("outputs/pheno.Rds")

pheno <- filter(pheno, !is.na(cage.pen))

pheno$feed <- ifelse(grepl("kontroll", pheno$feed),
                     "control",
                     "treatment")

write.csv(pheno[!grepl("^ftir_", colnames(pheno))],
          file = "tables/supplementary_data3_phenotypes.csv",
          quote = FALSE,
          row.names = FALSE)