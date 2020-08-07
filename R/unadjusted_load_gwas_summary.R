
## Summarise results of conditional body weight GWAS

library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(tidyr)
library(patchwork)
library(purrr)


source("R/gwas_helper_functions.R")

filenames_unadjusted <- system("ls gwas/*_unadjusted_load/output/*.assoc.txt",
                               intern = TRUE)
names(filenames_unadjusted) <- str_match(filenames_unadjusted, "all|cage|pen")[,1]

gwas_unadjusted <- map_dfr(filenames_unadjusted,
                           read_tsv,
                           col_types = "ccnnccnnnnnnnnn",
                           .id = "scan")

filenames_adjusted <- system("ls gwas/*load_N/output/*.assoc.txt",
                             intern = TRUE)
names(filenames_adjusted) <- str_match(filenames_adjusted, "all|cage|pen")[,1]

gwas_adjusted <- map_dfr(filenames_adjusted,
                         read_tsv,
                         col_types = "ccnnccnnnnnnnnn",
                         .id = "scan")

filenames_weight <- system("ls gwas/*weight/output/*.assoc.txt",
                             intern = TRUE)
names(filenames_weight) <- str_match(filenames_weight, "all|cage|pen")[,1]

gwas_weight <- map_dfr(filenames_weight,
                       read_tsv,
                       col_types = "ccnnccnnnnnnnnn",
                       .id = "scan")


## Set up chromosome lengths and marker positions for Manhattan plots

chr_lengths <- summarise(group_by(gwas_unadjusted, chr), length = unique(max(ps)))
preferred_order <- c(1:28, 30:31, 33, "Un_NW_020110160v1", "Un_NW_020110165v1", "Z")

chr_lengths <- chr_lengths[match(preferred_order, chr_lengths$chr),]

gwas_unadjusted$global_pos <- flatten_coordinates(gwas_unadjusted$chr,
                                                  gwas_unadjusted$ps,
                                                  chr_lengths)

gwas_adjusted$global_pos <- flatten_coordinates(gwas_adjusted$chr,
                                                gwas_adjusted$ps,
                                                chr_lengths)

gwas_weight$global_pos <- flatten_coordinates(gwas_weight$chr,
                                              gwas_weight$ps,
                                              chr_lengths)


## Find suggestive hits in default adjusted scan

suggestive_adjusted <- filter(gwas_adjusted, p_wald < 1e-4)

significant_weight <- filter(gwas_weight, p_wald < 5e-8)


chr_breaks <- summarise(group_by(gwas_unadjusted, chr),
                        global_pos_break = min(global_pos))

chr_breaks <- chr_breaks[match(preferred_order, chr_breaks$chr),]
chr_breaks$chr_masked <- chr_breaks$chr
chr_breaks$chr_masked[11:33] <- ""


formatting <- list(theme_bw(),
                   theme(legend.position = "none",
                         panel.grid = element_blank()),
                   scale_colour_manual(values = c("grey", "black")),
                   geom_hline(yintercept = 4, colour = "blue", type = 2),
                   ylab("-log(p)"))

plot_manhattan_cage <- plot_manhattan(filter(gwas_unadjusted, scan == "cage")) +
    geom_point(aes(x = global_pos,
                   y = 0),
               data = filter(suggestive_adjusted,
                             scan == "cage"),
               colour = "red",
               shape = 17) +
    geom_point(aes(x = global_pos,
                   y = 0),
               data = filter(significant_weight,
                             scan == "cage"),
               colour = "blue",
               shape = 17) +
    formatting +
    ggtitle("CAGE")
                         
plot_manhattan_pen <-  plot_manhattan(filter(gwas_unadjusted, scan == "pen")) +
                                        geom_point(aes(x = global_pos,
                                                       y = 0),
                                                   data = filter(suggestive_adjusted,
                                                                 scan == "pen"),
                                                   colour = "red",
                                                   shape = 17) +
                                        geom_point(aes(x = global_pos,
                                                       y = 0),
                                                   data = filter(significant_weight,
                                                                 scan == "pen"),
                                                   colour = "blue",
                                                   shape = 17) +
                                        formatting +
                                        ggtitle("PEN")

plot_manhattan_all <- plot_manhattan(filter(gwas_unadjusted, scan == "all")) +
    geom_point(aes(x = global_pos,
                   y = 0),
               data = filter(significant_weight,
                             scan == "all"),
               colour = "blue",
               shape = 17) +
    formatting +
    ggtitle("JOINT")



## Comparison between scans

comparison_all <- rbind(transform(filter(gwas_unadjusted, scan == "all"),
                                  case = "unadjusted"),
                        transform(filter(gwas_adjusted, scan == "all"),
                                  case = "adjusted"))


comparison_all_wide <- pivot_wider(comparison_all[, c("chr", "rs", "beta", "p_wald", "case")],
                                   values_from = c("beta", "p_wald"),
                                   names_from = "case")

formatting_comparison <- list(theme_bw(),
                              theme(panel.grid = element_blank()),
                              geom_abline(intercept = 0,
                                          slope = 1,
                                          colour = "blue"))


plot_comparison_beta <- qplot(comparison_all_wide$beta_adjusted,
                              comparison_all_wide$beta_unadjusted,
                              ylab = "Marker effect unadjusted",
                              xlab = "Marker effect adjusted") +
    formatting_comparison
                        
plot_comparison_p <- qplot(-log10(comparison_all_wide$p_wald_adjusted),
                           -log10(comparison_all_wide$p_wald_unadjusted),
                           ylab = "-log(p) unadjusted",
                           xlab = "-log(p) adjusted") +
    formatting_comparison



## Combined figures

plot_unadjusted <- ((plot_manhattan_cage /
                        plot_manhattan_pen /
                        plot_manhattan_all ) |
    (plot_comparison_beta /
         plot_comparison_p)) +
    plot_annotation(title = "Bone strength GWAS, not adjusted for body weight")



pdf("figures/plot_unadjusted_gwas.pdf")
print(plot_unadjusted)
dev.off()
