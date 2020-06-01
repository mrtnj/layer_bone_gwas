
## Analyse simulation results

library(dplyr)
library(GenomicRanges)
library(ggplot2)
library(readr)


source("R/gwas_helper_functions.R")

gwas <- read_tsv("simulations/output/simulation1.assoc.txt")

qtl <- read_delim("simulations/simulation1_qtl.txt", delim = " ")

qtl <- filter(qtl, frequency > 0)


chr_lengths <- data.frame(chr = 1:10,
                          length = rep(10^8, 10))

gwas$global_pos <- flatten_coordinates(gwas$chr,
                                       gwas$ps,
                                       chr_lengths)

chr_breaks <- summarise(group_by(gwas, chr),
                        global_pos_break = min(global_pos))

chr_breaks$chr_masked <- chr_breaks$chr



plot_manhattan <- plot_manhattan(gwas)



suggestive <- filter(gwas, p_wald < 10^-6)

suggestive_ranges <- GRanges(seqnames = suggestive$chr,
                             ranges = IRanges(suggestive$ps - 0.5e6,
                                              suggestive$ps + 0.5e6))

suggestive_regions <- reduce(suggestive_ranges)


qtl_ranges <- GRanges(seqnames = qtl$chr,
                      ranges = IRanges(qtl$pos, qtl$pos))
