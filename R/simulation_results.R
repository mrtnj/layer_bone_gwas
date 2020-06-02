
## Analyse simulation results

library(dplyr)
library(GenomicRanges)
library(ggplot2)
library(readr)


source("R/gwas_helper_functions.R")
source("R/simulation_helper_functions.R")

gwas_joint <- read_tsv("simulations/shared/output/shared1_joint.assoc.txt")
gwas_e1 <- read_tsv("simulations/shared/output/shared1_e1.assoc.txt")
gwas_e2 <- read_tsv("simulations/shared/output/shared1_e2.assoc.txt")

qtl <- read_delim("simulations/shared/shared1_qtl.txt", delim = " ")

qtl <- filter(qtl, frequency > 0)


chr_lengths <- data.frame(chr = 1:10,
                          length = rep(10^8, 10))

gwas$global_pos <- flatten_coordinates(gwas$chr,
                                       gwas$ps,
                                       chr_lengths)

chr_breaks <- summarise(group_by(gwas, chr),
                        global_pos_break = min(global_pos))

chr_breaks$chr_masked <- chr_breaks$chr



plot_joint <- plot_manhattan(gwas)



qtl$detected_joint <- find_detected_qtl(gwas_joint,
                                        threshold = 1e-6,
                                        flank = 0.5e6,
                                        qtl)

qtl$detected_e1 <- find_detected_qtl(gwas_e1,
                                     threshold = 1e-6,
                                     flank = 0.5e6,
                                     qtl)

qtl$detected_e2 <- find_detected_qtl(gwas_e2,
                                     threshold = 1e-6,
                                     flank = 0.5e6,
                                     qtl)

