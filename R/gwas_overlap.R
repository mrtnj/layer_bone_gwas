
## Read suggestive GWAS hits and make graphs

library(biomaRt)
library(dplyr)
library(GenomicRanges)
library(readr)
library(stringr)


gwas <- readRDS("outputs/gwas.Rds")


significant <- filter(gwas, p_wald < 5e-8)
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


bw_significant_regions <- reduce_region(filter(significant, scan_name_without_group == "weight"))



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



## Genes close to body weight regions

ensembl <- useMart("ensembl",
                   dataset = "ggallus_gene_ensembl")

get_region_genes <- function(regions) {

    getBM(attributes = c("ensembl_gene_id", "gene_biotype", "external_gene_name", "description"), 
          filters = "chromosomal_region",
          values = paste(seqnames(regions), ":",
                         start(regions), ":",
                         end(regions)),
          mart = ensembl)

}


## region_genes <- vector(mode = "list",
##                       length = length(bw_significant_regions))

##for (region_ix in 1:length(bw_significant_regions)) {

##    region_genes[[region_ix]] <- get_region_genes(bw_significant_regions[region_ix])

##}

##saveRDS(region_genes, "outputs/region_genes_bw.Rds")

region_genes <- readRDS("outputs/region_genes_bw.Rds")



## Animal QTLdb

qtl <- read_tsv("annotation/chickenqtldb2020-06_16.txt",
                col_names = FALSE,
                na = ".",
                comment = "#")

qtl$trait <- str_match(qtl$X9,
                       "Name=([^;]+);")[,2]
qtl$pmid <- str_match(qtl$X9,
                      "PUBMED_ID=([^;]+);")[,2]
qtl$markers <- str_match(qtl$X9,
                         "FlankMarker=([^;]+);")[,2]

qtl <- filter(qtl, grepl("Body weight", qtl$trait))

markers_split <- strsplit(qtl$markers, split = ",")

qtl$marker1 <- unlist(lapply(markers_split, "[", 1))
qtl$marker2 <- unlist(lapply(markers_split, "[", 2))

qtl$marker2 <- ifelse(is.na(qtl$marker2),
                      qtl$marker1,
                      qtl$marker2)

qtl$chr <- sub(qtl$X1,
               pattern = "Chr.",
               replacement = "")


variation <- read_tsv("annotation/gallus_gallus.gvf",
                      col_names = FALSE,
                      comment = "#")
variation <- filter(variation, X1 %in% c(4, 6, 27))

variation$id <- str_match(variation$X9, "Dbxref=dbSNP_150:([^;]+);")[,2]

variant_position <- variation[, c("id", "X4")]
colnames(variant_position) <- c("marker", "pos")


qtl <- inner_join(qtl, variant_position, by = c("marker1" = "marker"))
qtl <- inner_join(qtl, variant_position, by = c("marker2" = "marker"))


qtl_ranges <- GRanges(seqnames = qtl$chr,
                      ranges = IRanges(qtl$pos.x, qtl$pos.y))


qtldb_bw_overlap <- qtl[queryHits(findOverlaps(qtl_ranges, bw_significant_regions)),]
