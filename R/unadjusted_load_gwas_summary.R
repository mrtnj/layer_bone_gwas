
## Summarise results of conditional body weight GWAS

library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(purrr)


source("R/gwas_helper_functions.R")

filenames <- system("ls gwas/*_unadjusted_load/output/*.assoc.txt",
                    intern = TRUE)
names(filenames) <- str_match(filenames, "all|cage|pen")[,1]

gwas <- map_dfr(filenames,
                read_tsv,
                col_types = "ccnnccnnnnnnnnn",
                .id = "scan")


## Set up chromosome lengths and marker positions for Manhattan plots

chr_lengths <- summarise(group_by(gwas, chr), length = unique(max(ps)))
preferred_order <- c(1:28, 30:31, 33, "Un_NW_020110160v1", "Un_NW_020110165v1", "Z")

chr_lengths <- chr_lengths[match(preferred_order, chr_lengths$chr),]

gwas$global_pos <- flatten_coordinates(gwas$chr, gwas$ps, chr_lengths)


chr_breaks <- summarise(group_by(gwas, chr),
                        global_pos_break = min(global_pos))

chr_breaks <- chr_breaks[match(preferred_order, chr_breaks$chr),]
chr_breaks$chr_masked <- chr_breaks$chr
chr_breaks$chr_masked[11:33] <- ""


plot_manhattan(filter(gwas, scan == "all"))
plot_manhattan(filter(gwas, scan == "cage"))
plot_manhattan(filter(gwas, scan == "pen"))
