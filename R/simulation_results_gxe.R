
## Analyse simulation results

library(dplyr)
library(GenomicRanges)
library(ggplot2)
library(readr)


source("R/gwas_helper_functions.R")

gwas_joint <- read_tsv("simulations/gxe/output/gxe1_joint.assoc.txt")
gwas_e1 <- read_tsv("simulations/gxe/output/gxe1_e1.assoc.txt")
gwas_e2 <- read_tsv("simulations/gxe/output/gxe1_e2.assoc.txt")

qtl <- read_delim("simulations/gxe/gxe1_qtl.txt", delim = " ")

qtl <- filter(qtl, frequency > 0)


chr_lengths <- data.frame(chr = 1:10,
                          length = rep(10^8, 10))

gwas_joint$global_pos <- flatten_coordinates(gwas_joint$chr,
                                             gwas_joint$ps,
                                             chr_lengths)

gwas_e1$global_pos <- flatten_coordinates(gwas_e1$chr,
                                          gwas_e1$ps,
                                          chr_lengths)

gwas_e2$global_pos <- flatten_coordinates(gwas_e2$chr,
                                          gwas_e2$ps,
                                          chr_lengths)

chr_breaks <- summarise(group_by(gwas_joint, chr),
                        global_pos_break = min(global_pos))

chr_breaks$chr_masked <- chr_breaks$chr



plot_joint <- plot_manhattan(gwas_joint)
plot_e1 <- plot_manhattan(gwas_e1)
plot_e2 <- plot_manhattan(gwas_e2)


suggestive <- filter(gwas, p_wald < 10^-6)

suggestive_ranges <- GRanges(seqnames = suggestive$chr,
                             ranges = IRanges(suggestive$ps - 0.5e6,
                                              suggestive$ps + 0.5e6))

suggestive_regions <- reduce(suggestive_ranges)


qtl_ranges <- GRanges(seqnames = qtl$chr,
                      ranges = IRanges(qtl$pos, qtl$pos))
