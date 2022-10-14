
## Look at pre-defined candidate loci

library(dplyr)
library(readr)

candidate_loci <- read_csv("tables/candidate_regions.csv")

candidate_loci$start_expanded <- candidate_loci$start - 50e3
candidate_loci$end_expanded <- candidate_loci$end + 50e3



## GWAS resuls

gwas_pen_bovans_load <- readRDS("outputs/gwas_pen_bovans_load_coord.Rds")
gwas_cage_bovans_load <- readRDS("outputs/gwas_cage_bovans_load_coord.Rds")
gwas_all_bovans_load <- readRDS("outputs/gwas_all_bovans_load_coord.Rds")

gwas_pen_lsl_load <- readRDS("outputs/gwas_pen_lsl_load_coord.Rds")
gwas_cage_lsl_load <- readRDS("outputs/gwas_cage_lsl_load_coord.Rds")
gwas_all_lsl_load <- readRDS("outputs/gwas_all_lsl_load_coord.Rds")


gwas_markers <- unique(gwas_all_bovans_load[, c("chr", "marker_id", "ps")])

candidate_markers <- vector(mode = "list",
                            length = nrow(candidate_loci))

for (candidate_ix in 1:nrow(candidate_loci)) {

    candidate_markers[[candidate_ix]] <-
        gwas_markers[gwas_markers$chr == candidate_loci$chr[candidate_ix] &
                     gwas_markers$ps >= candidate_loci$start_expanded[candidate_ix] &
                     gwas_markers$ps <= candidate_loci$end_expanded[candidate_ix],]

}

candidates <- unique(Reduce(rbind, candidate_markers)$marker_id)


candidate_pen_bovans_load <- filter(gwas_pen_bovans_load,
                                    marker_id %in% candidates)

candidate_cage_bovans_load <- filter(gwas_cage_bovans_load,
                                     marker_id %in% candidates)

candidate_all_bovans_load <- filter(gwas_all_bovans_load,
                                    marker_id %in% candidates)


candidate_pen_lsl_load <- filter(gwas_pen_lsl_load,
                                 marker_id %in% candidates)

candidate_cage_lsl_load <- filter(gwas_cage_lsl_load,
                                  marker_id %in% candidates)

candidate_all_lsl_load <- filter(gwas_all_lsl_load,
                                 marker_id %in% candidates)


combined <- rbind(transform(candidate_cage_bovans_load,
                            scan_name = "Bone strength (Bovans, CAGE)"),
                  transform(candidate_pen_bovans_load,
                            scan_name = "Bone strength (Bovans, PEN)"),
                  transform(candidate_all_bovans_load,
                            scan_name = "Bone strength (Bovans, JOINT)"),
                  transform(candidate_cage_lsl_load,
                            scan_name = "Bone strength (LSL, CAGE)"),
                  transform(candidate_pen_lsl_load, 
                            scan_name = "Bone strength (LSL, PEN)"),
                  transform(candidate_all_lsl_load,
                            scan_name = "Bone strength (LSL, JOINT)"))



suggestive_candidates <- filter(combined, p < 1e-2)



write.csv(suggestive_candidates[order(as.numeric(suggestive_candidates$chr), suggestive_candidates$ps),],
          file = "tables/candidate_gwas.csv",
          quote = TRUE,
          row.names = FALSE)

