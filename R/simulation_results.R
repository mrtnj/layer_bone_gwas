
## Analyse simulation results

library(dplyr)
library(GenomicRanges)
library(ggplot2)
library(patchwork)
library(readr)


source("R/gwas_helper_functions.R")
source("R/simulation_helper_functions.R")


files_joint <- paste("simulations/shared/output/shared", 1:10,
                     "_joint.assoc.txt", sep = "")
files_e1 <- paste("simulations/shared/output/shared", 1:10,
                  "_e1.assoc.txt", sep = "")
files_e2 <- paste("simulations/shared/output/shared", 1:10,
                  "_e2.assoc.txt", sep = "")
files_qtl <- paste("simulations/shared/shared", 1:10,
                   "_qtl.txt", sep = "")


detected_qtl <- vector(length = 10,
                       mode = "list")

for (rep_ix in 1:10) {

    gwas_joint <- read_tsv(files_joint[rep_ix])
    gwas_e1 <- read_tsv(files_e1[rep_ix])
    gwas_e2 <- read_tsv(files_e2[rep_ix])
    
    qtl <- read_delim(files_qtl[rep_ix], delim = " ")
    
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
    
    pdf(paste("simulations/shared/shared", rep_ix, "_manhattan.pdf", sep = ""))
    print(plot_joint / plot_e1 / plot_e2)
    dev.off()
    
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
    
    detected_qtl[[rep_ix]] <- qtl
    
    write_tsv(qtl,
              paste("simulations/shared/shared", rep_ix, "_detected_qtl.txt", sep = ""))
    
    ## Scan comparison
    
    gwas_e1$scan_name <- "environment1"
    gwas_e2$scan_name <- "environment2"
    
    comparison <- pivot_wider(rbind(gwas_e1, gwas_e2)[, c("rs", "p_wald", "beta", "scan_name")],
                              values_from = c("beta", "p_wald"),
                              names_from = "scan_name")
    
    plot_comparison_p <- qplot(x = -log10(p_wald_environment1),
                               y = -log10(p_wald_environment2),
                               data = filter(comparison,
                                             p_wald_environment1 < 1e-3 |
                                                 p_wald_environment2 < 1e-3))
    
    plot_comparison_b <- qplot(x = beta_environment1,
                               y = beta_environment2,
                               data = filter(comparison,
                                             p_wald_environment1 < 1e-3 |
                                                 p_wald_environment2 < 1e-3)) +
        geom_hline(yintercept = 0, colour = "red", linetype = 2) +
        geom_vline(xintercept = 0, colour = "red", linetype = 2)
    
    pdf(paste("simulations/shared/shared", rep_ix, "_scan_comparison.pdf", sep = ""))
    print(plot_comparison_p / plot_comparison_b)
    dev.off()
}


detected_qtl_df <- Reduce(rbind, detected_qtl)
