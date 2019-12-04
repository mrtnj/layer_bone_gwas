
## Look at pre-defined candidate loci

library(dplyr)
library(readr)


candidate_geno <- read.table("mdb_dump/genotyping.txt",
                             header = TRUE,
                             stringsAsFactors = FALSE,
                             dec = ",",
                             sep = "\t")


write(colnames(candidate_geno)[-1], "outputs/candidate_markers.txt", sep = "\n")



## Load lifted positions from BioMart and manual lifting

mart <- read_tsv("candidate_regions/mart_export_candidate_markers.txt")
colnames(mart) <- c("id", "source", "chr", "start", "end")


extra <- read_csv("candidate_regions/candidate_regions_Galgal6.csv")


candidate_loci <- rbind(mart[,c("chr", "start", "end")],
                        extra)
candidate_loci$start_expanded <- candidate_loci$start - 50e3
candidate_loci$end_expanded <- candidate_loci$end + 50e3



## GWAS resuls

gwas <- readRDS("outputs/gwas.Rds")

gwas_markers <- unique(gwas[, c("chr", "rs", "ps")])

candidate_markers <- vector(mode = "list",
                            length = nrow(candidate_loci))

for (candidate_ix in 1:nrow(candidate_loci)) {

    candidate_markers[[candidate_ix]] <-
        gwas_markers[gwas_markers$chr == candidate_loci$chr[candidate_ix] &
                     gwas_markers$ps >= candidate_loci$start_expanded[candidate_ix] &
                     gwas_markers$ps <= candidate_loci$end_expanded[candidate_ix],]

}

candidates <- unique(Reduce(rbind, candidate_markers)$rs)


candidate_gwas <- filter(gwas, rs %in% candidates)



as.data.frame(filter(candidate_gwas, p_wald < 1e-2 & scan_name %in% c("pen_load_adj", "cage_load_adj", "all_load_adj")))
