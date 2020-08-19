
library(dplyr)
library(egg)
library(ggplot2)
library(patchwork)
library(readr)
library(stringr)
library(tidyr)

source("R/gwas_helper_functions.R")


gwas_pen_load <- readRDS("gwas/hglm_gwas_pen_load.Rds")
gwas_cage_load <- readRDS("gwas/hglm_gwas_cage_load.Rds")
gwas_all_load <- readRDS("gwas/hglm_gwas_all_load.Rds")

gwas_pen_weight <- readRDS("gwas/hglm_gwas_pen_weight.Rds")
gwas_cage_weight <- readRDS("gwas/hglm_gwas_cage_weight.Rds")
gwas_all_weight <- readRDS("gwas/hglm_gwas_all_weight.Rds")

gwas_all_weight_conditional <- readRDS("gwas/hglm_gwas_all_weight_conditional.Rds")


clean_marker_id <- function(marker_id) {
    sub(marker_id, pattern = "_[ATGC]", replacement = "")
}

gwas_pen_load$marker_id <- clean_marker_id(gwas_pen_load$marker_id)
gwas_cage_load$marker_id <- clean_marker_id(gwas_cage_load$marker_id)
gwas_all_load$marker_id <- clean_marker_id(gwas_all_load$marker_id) 

gwas_pen_weight$marker_id <- clean_marker_id(gwas_pen_weight$marker_id)
gwas_cage_weight$marker_id <- clean_marker_id(gwas_cage_weight$marker_id)
gwas_all_weight$marker_id <- clean_marker_id(gwas_all_weight$marker_id)

gwas_all_weight_conditional$marker_id <- clean_marker_id(gwas_all_weight_conditional$marker_id)

map <- read_delim("gwas/all.map",
                  delim = " ",
                  col_types = "ccnn",
                  col_names = FALSE)[,-3]

colnames(map) <- c("chr", "marker_id", "ps")

gwas_pen_load <- inner_join(map, gwas_pen_load)
gwas_cage_load <- inner_join(map, gwas_cage_load)
gwas_all_load <- inner_join(map, gwas_all_load)

gwas_pen_weight <- inner_join(map, gwas_pen_weight)
gwas_cage_weight <- inner_join(map, gwas_cage_weight)
gwas_all_weight <- inner_join(map, gwas_all_weight)

gwas_all_weight_conditional <- inner_join(map, gwas_all_weight_conditional)


## Set up chromosome lengths and marker positions for Manhattan plots

chr_lengths <- summarise(group_by(gwas_all_weight, chr), length = unique(max(ps)))
preferred_order <- c(1:28, 30:31, 33, "Z")

chr_lengths <- chr_lengths[match(preferred_order, chr_lengths$chr),]

gwas_pen_load$global_pos <- flatten_coordinates(gwas_pen_load$chr, gwas_pen_load$ps, chr_lengths)
gwas_cage_load$global_pos <- flatten_coordinates(gwas_cage_load$chr, gwas_cage_load$ps, chr_lengths)
gwas_all_load$global_pos <- flatten_coordinates(gwas_all_load$chr, gwas_all_load$ps, chr_lengths)

gwas_pen_weight$global_pos <- flatten_coordinates(gwas_pen_weight$chr, gwas_pen_weight$ps, chr_lengths)
gwas_cage_weight$global_pos <- flatten_coordinates(gwas_cage_weight$chr, gwas_cage_weight$ps, chr_lengths)
gwas_all_weight$global_pos <- flatten_coordinates(gwas_all_weight$chr, gwas_all_weight$ps, chr_lengths)

gwas_all_weight_conditional$global_pos <- flatten_coordinates(gwas_all_weight_conditional$chr,
                                                              gwas_all_weight_conditional$ps,
                                                              chr_lengths)

chr_breaks <- summarise(group_by(gwas_all_load, chr),
                        global_pos_break = min(global_pos))

chr_breaks <- chr_breaks[match(preferred_order, chr_breaks$chr),]
chr_breaks$chr_masked <- chr_breaks$chr
chr_breaks$chr_masked[11:32] <- ""


formatting <- list(geom_hline(yintercept = -log10(5e-8),
                              colour = "red",
                              linetype = 2),
                   geom_hline(yintercept = -log10(1e-4),
                              colour = "blue",
                              linetype = 2),
                   theme_bw(),
                   theme(panel.grid = element_blank(),
                         legend.position = "none"),
                   scale_colour_manual(values = c("grey", "black")),
                   xlab(""),
                   ylab(""))


## Main Manhattan plot of breaking strenght and body weight scans

plot_manhattan_cage_load <- plot_manhattan(gwas_cage_load,
                                          p = "p") +
    formatting +
    ylim(0, 10) +
    ggtitle("Bone breaking strength (CAGE)")

plot_manhattan_pen_load <- plot_manhattan(gwas_pen_load,
                                           p = "p") +
    formatting +
    ylim(0, 10) +
    ggtitle("Bone breaking strength (PEN)")

plot_manhattan_all_load <- plot_manhattan(gwas_all_load,
                                          p = "p") +
    formatting +
    ylim(0, 10) +
    ggtitle("Bone breaking strength (JOINT)")



plot_manhattan_cage_weight <- plot_manhattan(gwas_cage_weight,
                                             p = "p") +
    formatting +
    ylim(0, 17) +
    ggtitle("Body weight (CAGE)")

plot_manhattan_pen_weight <- plot_manhattan(gwas_pen_weight,
                                          p = "p") +
    formatting +
    ylim(0, 17) +
    ggtitle("Body weight (PEN)")

plot_manhattan_all_weight <- plot_manhattan(gwas_all_weight,
                                          p = "p") +
    formatting +
    ylim(0, 17) +
    ggtitle("Body weight (JOINT)")


plot_manhattan_combined <- ggarrange(plot_manhattan_cage_load,
                                     plot_manhattan_pen_load,
                                     plot_manhattan_all_load,
                                     plot_manhattan_cage_weight,
                                     plot_manhattan_pen_weight,
                                     plot_manhattan_all_weight,
                                     left = "Negative logarithm of p-value",
                                     ncol = 2,
                                     byrow = FALSE)


pdf("figures/plot_manhattan_hglm.pdf")
print(plot_manhattan_combined)
dev.off()


plot_qq_cage_load <- plot_qq(gwas_cage_load$p) +
    ggtitle("Bone breaking strength (CAGE)")

plot_qq_pen_load <- plot_qq(gwas_pen_load$p) +
    ggtitle("Bone breaking strength (PEN)")

plot_qq_all_load <- plot_qq(gwas_all_load$p) +
    ggtitle("Bone breaking strength (JOINT)")



plot_qq_cage_weight <- plot_qq(gwas_cage_weight$p) +
    ggtitle("Body weight (CAGE)")

plot_qq_pen_weight <- plot_qq(gwas_pen_weight$p) +
    ggtitle("Body weight (PEN)")

plot_qq_all_weight <- plot_qq(gwas_all_weight$p) +
    ggtitle("Body weight (JOINT)")



plot_qq_combined <- ggarrange(plot_qq_cage_load,
                              plot_qq_pen_load,
                              plot_qq_all_load,
                              plot_qq_cage_weight,
                              plot_qq_pen_weight,
                              plot_qq_all_weight,
                              ncol = 2,
                              byrow = FALSE)

pdf("figures/plot_qq_hglm.pdf")
print(plot_qq_combined)
dev.off()


## Conditional scan

combined_conditional <- rbind(transform(filter(gwas_all_weight, chr == 4),
                                        scan = "Default"),
                              transform(filter(gwas_all_weight_conditional, chr == 4),
                                        scan = "Conditional"))

plot_conditional <- qplot(x = ps/1e6, y = -log10(p),
                          colour = scan,
                          data = combined_conditional) +
    formatting +
    ggtitle("Conditional GWAS of body weight on chromosome 4")


pdf("figures/plot_qq_hglm.pdf")
print(plot_qq_combined)
dev.off()


## Comparison of systems

gwas_load_systems <- rbind(transform(gwas_pen_load,
                                     scan_name = "pen"),
                           transform(gwas_cage_load,
                                     scan_name = "cage"))

load_comparison <- pivot_wider(gwas_load_systems[, c("marker_id", "ps", "p", "estimates", "scan_name")],
                               values_from = c("estimates", "p"),
                               names_from = "scan_name")


gwas_weight_systems <- rbind(transform(gwas_pen_weight,
                                       scan_name = "pen"),
                             transform(gwas_cage_weight,
                                       scan_name = "cage"))

weight_comparison <- pivot_wider(gwas_weight_systems[, c("marker_id", "ps", "p", "estimates", "scan_name")],
                               values_from = c("estimates", "p"),
                               names_from = "scan_name")

formatting_comparison <- list(theme_bw(),
                              theme(panel.grid = element_blank()),
                              xlab("CAGE"),
                              ylab("PEN"))


plot_load_comparison_p <- qplot(x = -log10(p_cage),
                                y = -log10(p_pen),
                                data = filter(load_comparison,
                                              p_cage < 1e-3 |
                                                  p_pen < 1e-3)) +
    ggtitle("Bone strength \nNegative logarithm of p-value") +
    formatting_comparison +
    xlim(0, 10) +
    ylim(0, 10)

plot_load_comparison_beta <- qplot(x = estimates_cage,
                                   y = estimates_pen,
                                   data = filter(load_comparison,
                                                 p_cage < 1e-3 |
                                                     p_pen < 1e-3)) +
    geom_hline(yintercept = 0, colour = "red", linetype = 2) +
    geom_vline(xintercept = 0, colour = "red", linetype = 2) +
    ggtitle("Bone strength \nEstimated marker effect") +
    formatting_comparison


plot_weight_comparison_p <- qplot(x = -log10(p_cage),
                                  y = -log10(p_pen),
                                  data = filter(weight_comparison,
                                                p_cage < 1e-3 |
                                                    p_pen < 1e-3)) +
    ggtitle("Body weight \nNegative logarithm of p-value") +
    formatting_comparison +
    xlim(0, 10) +
    ylim(0, 10)

plot_weight_comparison_beta <- qplot(x = estimates_cage,
                                     y = estimates_pen,
                                     data = filter(weight_comparison,
                                                   p_cage < 1e-3 |
                                                       p_pen < 1e-3)) +
    geom_hline(yintercept = 0, colour = "red", linetype = 2) +
    geom_vline(xintercept = 0, colour = "red", linetype = 2) +
    ggtitle("Body weight \nEstimated marker effect") +
    formatting_comparison



plot_comparisons <- ggarrange(plot_load_comparison_p,
                              plot_load_comparison_beta,
                              plot_weight_comparison_p,
                              plot_weight_comparison_beta)



pdf("figures/plot_gwas_comparison_hglm.pdf")
print(plot_comparisons)
dev.off()


