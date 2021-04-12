
## Look at pre-defined candidate loci

library(dplyr)
library(readr)


## Candidate genotypes

candidate_geno <- read.table("data/genotyping.txt",
                             header = TRUE,
                             stringsAsFactors = FALSE,
                             dec = ",",
                             sep = "\t")


write(colnames(candidate_geno)[-1], "outputs/candidate_markers.txt", sep = "\n")


## Regions defined on Galgal4

regions <- read_csv("candidate_regions/candidate_regions_Galgal4.csv")

## Write output for coordinate conversion

write.table(data.frame(paste("chr", regions$chr, sep = ""),
                       regions$start - 1,
                       regions$start),
            file = "outputs/regions_Galgal4_start.bed",
            sep = "\t",
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE)

write.table(data.frame(paste("chr", regions$chr, sep = ""),
                       regions$end - 1,
                       regions$end),
            file = "outputs/regions_Galgal4_end.bed",
            sep = "\t",
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE)



## Load lifted positions from BioMart and manual lifting

bed_start <- read_tsv("candidate_regions/candidate_regions_LiftOver_start.bed",
                      col_names = FALSE)
bed_end <- read_tsv("candidate_regions/candidate_regions_LiftOver_end.bed",
                    col_names = FALSE)

candidate_loci <- data.frame(chr = sub(bed_start$X1, pattern = "chr", replacement = ""),
                             start = bed_start$X3,
                             end = bed_end$X3,
                             stringsAsFactors = FALSE)

extra <- read_csv("candidate_regions/candidate_regions_Galgal6.csv")

candidate_loci <- rbind(candidate_loci,
                        extra)


## Flip chr16
chr16_starts <- candidate_loci$start[candidate_loci$chr == 16]
candidate_loci$start[candidate_loci$chr == 16] <- candidate_loci$end[candidate_loci$chr == 16]
candidate_loci$end[candidate_loci$chr == 16] <- chr16_starts

candidate_loci$start_expanded <- candidate_loci$start - 50e3
candidate_loci$end_expanded <- candidate_loci$end + 50e3



## GWAS resuls

gwas_pen_load <- readRDS("outputs/gwas_pen_load_coord.Rds")
gwas_cage_load <- readRDS("outputs/gwas_cage_load_coord.Rds")
gwas_all_load <- readRDS("outputs/gwas_all_load_coord.Rds")

gwas_markers <- unique(gwas_all_load[, c("chr", "marker_id", "ps")])

candidate_markers <- vector(mode = "list",
                            length = nrow(candidate_loci))

for (candidate_ix in 1:nrow(candidate_loci)) {

    candidate_markers[[candidate_ix]] <-
        gwas_markers[gwas_markers$chr == candidate_loci$chr[candidate_ix] &
                     gwas_markers$ps >= candidate_loci$start_expanded[candidate_ix] &
                     gwas_markers$ps <= candidate_loci$end_expanded[candidate_ix],]

}

candidates <- unique(Reduce(rbind, candidate_markers)$marker_id)


candidate_pen_load <- filter(gwas_pen_load,
                             marker_id %in% candidates)

candidate_cage_load <- filter(gwas_cage_load,
                              marker_id %in% candidates)

candidate_all_load <- filter(gwas_all_load,
                             marker_id %in% candidates)



suggestive_candidates <- rbind(transform(filter(candidate_cage_load, p < 1e-2),
                                         scan_name = "Bone strength (CAGE)"),
                               transform(filter(candidate_pen_load, p < 1e-2),
                                         scan_name = "Bone strength (PEN)"),
                               transform(filter(candidate_all_load, p < 1e-2),
                                         scan_name = "Bone strength (JOINT)"))


write.csv(suggestive_candidates[order(as.numeric(suggestive_candidates$chr), suggestive_candidates$ps),],
          file = "tables/candidate_gwas.csv",
          quote = FALSE,
          row.names = FALSE)

write.csv(candidate_loci[order(as.numeric(candidate_loci$chr), candidate_loci$start), 1:3],
          file = "tables/candidate_regions.csv",
          quote = FALSE,
          row.names = FALSE)
