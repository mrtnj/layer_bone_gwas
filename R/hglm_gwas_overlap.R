
## Find regions around GWAS hits for summary table

library(biomaRt)
library(dplyr)
library(GenomicRanges)
library(readr)
library(stringr)
library(purrr)

## Bone strength

gwas_pen_load <- readRDS("outputs/gwas_pen_load_coord.Rds")
gwas_pen_load$pretty_name <- "bone strength"
gwas_pen_load$scan <- "pen"

gwas_cage_load <- readRDS("outputs/gwas_cage_load_coord.Rds")
gwas_cage_load$pretty_name <- "bone strength"
gwas_cage_load$scan <- "cage"

gwas_all_load <- readRDS("outputs/gwas_all_load_coord.Rds")
gwas_all_load$pretty_name <- "bone strength"
gwas_all_load$scan <- "all"


## Body weight

gwas_pen_weight <- readRDS("outputs/gwas_pen_weight_coord.Rds")
gwas_pen_weight$pretty_name <- "body weight"
gwas_pen_weight$scan <- "pen"

gwas_cage_weight <- readRDS("outputs/gwas_cage_weight_coord.Rds")
gwas_cage_weight$pretty_name <- "body weight"
gwas_cage_weight$scan <- "cage"

gwas_all_weight <- readRDS("outputs/gwas_all_weight_coord.Rds")
gwas_all_weight$pretty_name <- "body weight"
gwas_all_weight$scan <- "all"


## pQCT & TGA

gwas_pen <- readRDS("outputs/gwas_pen_bone_phenotypes_coord.Rds")
gwas_cage <- readRDS("outputs/gwas_cage_bone_phenotypes_coord.Rds")
gwas_all <- readRDS("outputs/gwas_all_bone_phenotypes_coord.Rds")

gwas_pen$scan <- "pen"
gwas_cage$scan <- "cage"
gwas_all$scan <- "all"


cols <- c("pretty_name", "chr", "ps", "p", "scan")

combined <- rbind(gwas_pen_load[, cols],
                  gwas_cage_load[, cols],
                  gwas_all_load[, cols],
                  gwas_pen_weight[, cols],
                  gwas_cage_weight[, cols],
                  gwas_all_weight[, cols],
                  gwas_pen[, cols],
                  gwas_cage[, cols],
                  gwas_all[, cols])


significant <- filter(combined, p < 5e-8)

suggestive <- filter(combined, p < 1e-4)

reduce_region <- function(trait_markers) {
 
    marker_ranges <- GRanges(seqnames = trait_markers$chr,
                             ranges = IRanges(trait_markers$ps - 5e6,
                                              trait_markers$ps + 5e6))
    GenomicRanges::reduce(marker_ranges)
}


hit_traits <- unique(c(significant$pretty_name, suggestive$pretty_name))

significant_split <- lapply(hit_traits, function(trait) significant[significant$pretty_name == trait,])
suggestive_split <- lapply(hit_traits, function(trait) suggestive[suggestive$pretty_name == trait,])

names(significant_split) <- hit_traits
names(suggestive_split) <- hit_traits

significant_regions <- lapply(significant_split, reduce_region)
suggestive_regions <- lapply(suggestive_split, reduce_region)

suggestive_regions_pruned <- mapply(function(ranges1, ranges2) ranges1[!ranges1 %over% ranges2],
                                    suggestive_regions,
                                    significant_regions,
                                    SIMPLIFY = FALSE)

map_dfr(significant_regions, as.data.frame, .id = "trait")

map_dfr(suggestive_regions, as.data.frame, .id = "trait")
