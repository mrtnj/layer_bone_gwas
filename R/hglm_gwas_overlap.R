
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
gwas_pen_load$scan <- "(PEN)"

gwas_cage_load <- readRDS("outputs/gwas_cage_load_coord.Rds")
gwas_cage_load$pretty_name <- "bone strength"
gwas_cage_load$scan <- "(CAGE)"

gwas_all_load <- readRDS("outputs/gwas_all_load_coord.Rds")
gwas_all_load$pretty_name <- "bone strength"
gwas_all_load$scan <- "(JOINT)"


## Body weight

gwas_pen_weight <- readRDS("outputs/gwas_pen_weight_coord.Rds")
gwas_pen_weight$pretty_name <- "body weight"
gwas_pen_weight$scan <- "(PEN)"

gwas_cage_weight <- readRDS("outputs/gwas_cage_weight_coord.Rds")
gwas_cage_weight$pretty_name <- "body weight"
gwas_cage_weight$scan <- "(CAGE)"

gwas_all_weight <- readRDS("outputs/gwas_all_weight_coord.Rds")
gwas_all_weight$pretty_name <- "body weight"
gwas_all_weight$scan <- "(JOINT)"


## pQCT & TGA

gwas_pen <- readRDS("outputs/gwas_pen_bone_phenotypes_coord.Rds")
gwas_cage <- readRDS("outputs/gwas_cage_bone_phenotypes_coord.Rds")
gwas_all <- readRDS("outputs/gwas_all_bone_phenotypes_coord.Rds")

gwas_pen$scan <- "(PEN)"
gwas_cage$scan <- "(CAGE)"
gwas_all$scan <- "(JOINT)"


cols <- c("pretty_name", "chr", "ps", "marker_id", "estimates", "LRT", "p", "scan")

combined <- rbind(gwas_pen_load[, cols],
                  gwas_cage_load[, cols],
                  gwas_all_load[, cols],
                  gwas_pen_weight[, cols],
                  gwas_cage_weight[, cols],
                  gwas_all_weight[, cols],
                  gwas_pen[, cols],
                  gwas_cage[, cols],
                  gwas_all[, cols])

combined$pretty_name <- paste(combined$pretty_name,
                              combined$scan)

write.csv(combined,
          "tables/supplementary_data1_gwas_summary_stats.csv",
          quote = FALSE,
          row.names = FALSE)


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

## We now have the regions that are significant/suggestive for each trait. Now
## we need to go back and find the underlying SNPs


## Function to extract overlapping SNPs from a GenomicRanges object of QTL regions from a scan

get_snps <- function(regions, scan_name, combined) {
    
    snps_for_this_trait <- combined[combined$pretty_name == scan_name,]
    
    snp_ranges <- GRanges(seqnames = snps_for_this_trait$chr,
                          ranges = IRanges(snps_for_this_trait$ps, snps_for_this_trait$ps),
                          mcols = snps_for_this_trait)
    
    n_qtl <- length(regions)
    
    qtl_ranges_list <- vector(mode = "list",
                              length = n_qtl)
    
    if (n_qtl > 0) {
        
        for (qtl_ix in 1:n_qtl) {
            qtl_ranges_list[[qtl_ix]] <- subsetByOverlaps(snp_ranges,
                                                          regions[qtl_ix])
            
        }
    }
    
    qtl_ranges_list
}


## Find the most significant SNP in each region, extract statistics

get_top_snps <- function(regions) {

    n_qtl <- length(regions)
    
    top_snps <- vector(mode = "list",
                       length = n_qtl)
    
    if (n_qtl > 0) {
        
        for (qtl_ix in 1:n_qtl) {
            region_df <- as.data.frame(regions[[qtl_ix]])
            top_snps[[qtl_ix]] <- region_df[region_df$mcols.p == min(region_df$mcols.p),]
            
            ## If there are more than one, pick the first
            if (nrow(top_snps[[qtl_ix]]) > 1) {
                top_snps[[qtl_ix]] <- top_snps[[qtl_ix]][1,]
            }
        }
    }
    
    top_snps
}


## Apply this function to all the significant/suggestive scans

snps_significant_regions <- vector(mode = "list",
                                   length = length(significant_regions))

top_snps_significant_regions <- vector(mode = "list",
                                       length = length(significant_regions))

for (scan_ix in 1:length(significant_regions)) {
    snps_significant_regions[[scan_ix]] <- get_snps(significant_regions[[scan_ix]],
                                                    names(significant_regions)[scan_ix],
                                                    combined)
    top_snps_significant_regions[[scan_ix]] <- get_top_snps(snps_significant_regions[[scan_ix]])
}


snps_suggestive_regions <- vector(mode = "list",
                                  length = length(suggestive_regions_pruned))

top_snps_suggestive_regions <- vector(mode = "list",
                                      length = length(suggestive_regions_pruned))

for (scan_ix in 1:length(suggestive_regions_pruned)) {
    snps_suggestive_regions[[scan_ix]] <- get_snps(suggestive_regions_pruned[[scan_ix]],
                                                   names(suggestive_regions_pruned)[scan_ix],
                                                   combined)
    top_snps_suggestive_regions[[scan_ix]] <- get_top_snps(snps_suggestive_regions[[scan_ix]])
}




## Make tables with top SNPs

make_table <- function(scans_list,
                       top_snps_list) {

    significant_regions_df <- vector(mode = "list",
                                     length = length(scans_list))
    
    for (scan_ix in 1:length(scans_list)) {
        
        regions <- as.data.frame(scans_list[[scan_ix]])
        regions$seqnames <- as.character(regions$seqnames)
        
        n_qtl <- nrow(regions)
        
        if (n_qtl > 0) {
            regions$trait <- names(scans_list)[scan_ix]
            top_snps <- top_snps_list[[scan_ix]]
            
            regions$lead_p <- 0
            regions$lead_pos <- 0
            
            for (region_ix in 1:n_qtl) {
                regions$lead_p[region_ix] <- top_snps[[region_ix]]$mcols.p
                regions$lead_pos[region_ix] <- top_snps[[region_ix]]$mcols.ps
            }
        }
        
        significant_regions_df[[scan_ix]] <- regions
    }
 
    Reduce(rbind, significant_regions_df)
}

significant_table <- make_table(significant_regions, top_snps_significant_regions)

suggestive_table <- make_table(suggestive_regions_pruned, top_snps_suggestive_regions)



significant_table <- significant_table[order(as.numeric(significant_table$seqnames),
                                             significant_table$lead_pos),]

suggestive_table <- suggestive_table[order(as.numeric(suggestive_table$seqnames),
                                     suggestive_table$lead_pos),]



cols <- c("trait", "seqnames", "lead_pos", "lead_p")

write.csv(significant_table[, cols],
          file = "tables/significant_gwas_table.csv",
          quote = TRUE,
          row.names = FALSE)

write.csv(suggestive_table[, cols],
          file = "tables/suggestive_gwas_table.csv",
          quote = TRUE,
          row.names = FALSE)

