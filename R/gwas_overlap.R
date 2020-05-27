
## Read suggestive GWAS hits and make graphs


library(dplyr)
library(GenomicRanges)
library(readr)


gwas <- readRDS("outputs/gwas.Rds")


suggestive <- filter(gwas, p_wald < 1e-4)


marker_ranges <- GRanges(seqnames = suggestive$chr,
                         ranges = IRanges(suggestive$ps - 0.5e6,
                                          suggestive$ps + 0.5e6),
                         mcols = as.data.frame(suggestive))      


reduce_region <- function(trait_markers) {
 
    marker_ranges <- GRanges(seqnames = trait_markers$chr,
                             ranges = IRanges(trait_markers$ps - 0.5e6,
                                              trait_markers$ps + 0.5e6))      
    reduce(marker_ranges)
}



## Pen load overlaps -- are there any?

pen_regions <- reduce_region(filter(suggestive, scan_name == "pen_load_N"))

print(subsetByOverlaps(marker_ranges, pen_regions))
## No



## Cage load overlaps -- are there any?

cage_regions <- reduce_region(filter(suggestive, scan_name == "cage_load_N"))

print(subsetByOverlaps(marker_ranges, cage_regions))
## No. Only a composition overlap from the pen


## Body weight overlaps


bw_regions <- reduce_region(filter(suggestive, scan_name_without_group == "weight"))

print(as.data.frame(subsetByOverlaps(marker_ranges, bw_regions)))

## Not much




## Table of all regions

suggestive_reduced <- do(group_by(suggestive, scan_name),
                         as.data.frame(reduce_region(.)))


suggestive_reduced_pretty_names <- inner_join(suggestive_reduced[, 1:4],
                                              unique(suggestive[, c("scan_name", "scan_name_pretty")]))


colnames(suggestive_reduced_pretty_names)[c(2, 5)] <- c("chromosome", "trait")

suggestive_reduced_pretty_names <-
    suggestive_reduced_pretty_names[order(as.numeric(suggestive_reduced_pretty_names$chromosome),
                                          suggestive_reduced_pretty_names$start,
                                          suggestive_reduced_pretty_names$trait),]


write_csv(suggestive_reduced_pretty_names[, c(5, 2:4)],
          "tables/supplementary_table_all_regions.csv")
