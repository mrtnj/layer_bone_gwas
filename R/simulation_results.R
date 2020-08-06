
## Analyse simulation results

library(dplyr)
library(GenomicRanges)
library(ggplot2)
library(patchwork)
library(readr)
library(tidyr)


source("R/gwas_helper_functions.R")
source("R/simulation_helper_functions.R")

n_rep <- 50

files_joint_shared <- paste("simulations/shared/output/shared", 1:n_rep,
                            "_joint.assoc.txt", sep = "")
files_e1_shared <- paste("simulations/shared/output/shared", 1:n_rep,
                         "_e1.assoc.txt", sep = "")
files_e2_shared <- paste("simulations/shared/output/shared", 1:n_rep,
                         "_e2.assoc.txt", sep = "")
files_qtl_shared <- paste("simulations/shared/shared", 1:n_rep,
                          "_qtl.txt", sep = "")

files_joint_gxe <- paste("simulations/gxe/output/gxe", 1:n_rep,
                         "_joint.assoc.txt", sep = "")
files_e1_gxe <- paste("simulations/gxe/output/gxe", 1:n_rep,
                      "_e1.assoc.txt", sep = "")
files_e2_gxe <- paste("simulations/gxe/output/gxe", 1:n_rep,
                      "_e2.assoc.txt", sep = "")
files_qtl_gxe <- paste("simulations/gxe/gxe", 1:n_rep,
                       "_qtl.txt", sep = "")


summarise_results <- function(files_joint,
                              files_e1,
                              files_e2,
                              files_qtl) {

    n_rep <- length(files_joint)
    
    detected_qtl <- vector(length = n_rep,
                           mode = "list")
    
    manhattan_plots <- vector(length = n_rep,
                              mode = "list")
    
    comparison_plots <- vector(length = n_rep,
                               mode = "list")
    
    overlapping_markers <- numeric(n_rep)
    cor_beta <- numeric(n_rep)
    
    for (rep_ix in 1:n_rep) {
        
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
        
        
        manhattan_plots[[rep_ix]] <- (plot_joint / plot_e1 / plot_e2)
        
        
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
        
        
        
        ## Scan comparison
        
        gwas_e1$scan_name <- "environment1"
        gwas_e2$scan_name <- "environment2"
        
        comparison <- pivot_wider(rbind(gwas_e1, gwas_e2)[, c("rs", "p_wald", "beta", "scan_name")],
                                  values_from = c("beta", "p_wald"),
                                  names_from = "scan_name")
        
        comparison_suggestive <- filter(comparison,
                                        p_wald_environment1 < 1e-3 |
                                            p_wald_environment2 < 1e-3)
        
        plot_comparison_p <- qplot(x = -log10(p_wald_environment1),
                                   y = -log10(p_wald_environment2),
                                   data = comparison_suggestive)
        
        plot_comparison_b <- qplot(x = beta_environment1,
                                   y = beta_environment2,
                                   data = comparison_suggestive) +
            geom_hline(yintercept = 0, colour = "red", linetype = 2) +
            geom_vline(xintercept = 0, colour = "red", linetype = 2)
        
        
        comparison_plots[[rep_ix]] <- (plot_comparison_p / plot_comparison_b)
        
        
        ## Number of shared markers and correlation b/w effects
        
        overlapping_markers[rep_ix] <- sum(comparison$p_wald_environment1 < 1e-3 &
                                               comparison$p_wald_environment2 < 1e-3,
                                           na.rm = TRUE)
        cor_beta[rep_ix] <- cor(comparison_suggestive$beta_environment1,
                                comparison_suggestive$beta_environment2,
                                use = "p")
    }
    
    
    detected_qtl_df <- Reduce(rbind, detected_qtl)
 
    list(detected_qtl = detected_qtl_df,
         manhattan_plots = manhattan_plots,
         comparison_plots = comparison_plots,
         overlapping_markers = overlapping_markers,
         cor_beta = cor_beta)   
    
}

results_shared <- summarise_results(files_joint_shared,
                                    files_e1_shared,
                                    files_e2_shared,
                                    files_qtl_shared)

results_gxe <- summarise_results(files_joint_gxe,
                                 files_e1_gxe,
                                 files_e2_gxe,
                                 files_qtl_gxe)

shared_power1 <- sum(results_shared$detected_qtl$percent_variance_explained > 1 &
                         results_shared$detected_qtl$detected_joint) /
    sum(results_shared$detected_qtl$percent_variance_explained > 1)

shared_power1_e1 <- sum(results_shared$detected_qtl$percent_variance_explained > 1 &
                         results_shared$detected_qtl$detected_e1) /
    sum(results_shared$detected_qtl$percent_variance_explained > 1)

gxe_power1_e1 <- sum(results_gxe$detected_qtl$percent_variance_explained1 > 1 &
                       results_gxe$detected_qtl$detected_e1) /
    sum(results_gxe$detected_qtl$percent_variance_explained1 > 1)

gxe_power1_e2 <- sum(results_gxe$detected_qtl$percent_variance_explained2 > 1 &
                         results_gxe$detected_qtl$detected_e2) /
    sum(results_gxe$detected_qtl$percent_variance_explained2 > 1)



plot_power_shared <- qplot(x = percent_variance_explained,
                           fill = detected_joint,
                           data = results_shared$detected_qtl) +
    annotate("text",
             label = paste("Power for loci > 1%: ", signif(shared_power1, 2)),
             x = 20, y = 175)



plot_power_gxe1 <- qplot(x = percent_variance_explained1,
                         fill = detected_e1,
                         data = results_gxe$detected_qtl) +
    annotate("text",
             label = paste("Power for loci > 1%: ", signif(gxe_power1_e1, 2)),
             x = 20, y = 175)

plot_power_gxe2 <- qplot(x = percent_variance_explained2,
                         fill = detected_e2,
                         data = results_gxe$detected_qtl) +
    annotate("text",
             label = paste("Power for loci > 1%: ", signif(gxe_power1_e2, 2)),
             x = 20, y = 175)

