
## Summarise the results of separate cross gwas

library(dplyr)
library(egg)
library(ggplot2)
library(patchwork)
library(readr)
library(stringr)
library(tidyr)
library(purrr)
library(patchwork)

source("R/hglm_helper_functions.R")
source("R/gwas_helper_functions.R")


gwas_pen_bovans_load <- readRDS("gwas/hglm_gwas_pen_bovans_load.Rds")
gwas_cage_bovans_load <- readRDS("gwas/hglm_gwas_cage_bovans_load.Rds")
gwas_all_bovans_load <- readRDS("gwas/hglm_gwas_all_bovans_load.Rds")

gwas_pen_lsl_load <- readRDS("gwas/hglm_gwas_pen_lsl_load.Rds")
gwas_cage_lsl_load <- readRDS("gwas/hglm_gwas_cage_lsl_load.Rds")
gwas_all_lsl_load <- readRDS("gwas/hglm_gwas_all_lsl_load.Rds")


gwas_pen_bovans_weight <- readRDS("gwas/hglm_gwas_pen_bovans_weight.Rds")
gwas_cage_bovans_weight <- readRDS("gwas/hglm_gwas_cage_bovans_weight.Rds")
gwas_all_bovans_weight <- readRDS("gwas/hglm_gwas_all_bovans_weight.Rds")

gwas_pen_lsl_weight <- readRDS("gwas/hglm_gwas_pen_lsl_weight.Rds")
gwas_cage_lsl_weight <- readRDS("gwas/hglm_gwas_cage_lsl_weight.Rds")
gwas_all_lsl_weight <- readRDS("gwas/hglm_gwas_all_lsl_weight.Rds")




clean_marker_id <- function(marker_id) {
    sub(marker_id, pattern = "_[ATGC]", replacement = "")
}

gwas_pen_bovans_load$marker_id <- clean_marker_id(gwas_pen_bovans_load$marker_id)
gwas_cage_bovans_load$marker_id <- clean_marker_id(gwas_cage_bovans_load$marker_id)
gwas_all_bovans_load$marker_id <- clean_marker_id(gwas_all_bovans_load$marker_id)

gwas_pen_lsl_load$marker_id <- clean_marker_id(gwas_pen_lsl_load$marker_id)
gwas_cage_lsl_load$marker_id <- clean_marker_id(gwas_cage_lsl_load$marker_id)
gwas_all_lsl_load$marker_id <- clean_marker_id(gwas_all_lsl_load$marker_id)


gwas_pen_bovans_weight$marker_id <- clean_marker_id(gwas_pen_bovans_weight$marker_id)
gwas_cage_bovans_weight$marker_id <- clean_marker_id(gwas_cage_bovans_weight$marker_id)
gwas_all_bovans_weight$marker_id <- clean_marker_id(gwas_all_bovans_weight$marker_id)

gwas_pen_lsl_weight$marker_id <- clean_marker_id(gwas_pen_lsl_weight$marker_id)
gwas_cage_lsl_weight$marker_id <- clean_marker_id(gwas_cage_lsl_weight$marker_id)
gwas_all_lsl_weight$marker_id <- clean_marker_id(gwas_all_lsl_weight$marker_id)




map <- read_delim("gwas/all.map",
                  delim = " ",
                  col_types = "ccnn",
                  col_names = FALSE)[,-3]

colnames(map) <- c("chr", "marker_id", "ps")

gwas_pen_bovans_load <- inner_join(map, gwas_pen_bovans_load)
gwas_cage_bovans_load <- inner_join(map, gwas_cage_bovans_load)
gwas_all_bovans_load <- inner_join(map, gwas_all_bovans_load)

gwas_pen_lsl_load <- inner_join(map, gwas_pen_lsl_load)
gwas_cage_lsl_load <- inner_join(map, gwas_cage_lsl_load)
gwas_all_lsl_load <- inner_join(map, gwas_all_lsl_load)


gwas_pen_bovans_weight <- inner_join(map, gwas_pen_bovans_weight)
gwas_cage_bovans_weight <- inner_join(map, gwas_cage_bovans_weight)
gwas_all_bovans_weight <- inner_join(map, gwas_all_bovans_weight)

gwas_pen_lsl_weight <- inner_join(map, gwas_pen_lsl_weight)
gwas_cage_lsl_weight <- inner_join(map, gwas_cage_lsl_weight)
gwas_all_lsl_weight <- inner_join(map, gwas_all_lsl_weight)


saveRDS(gwas_pen_bovans_load,
        file = "outputs/gwas_pen_bovans_load_coord.Rds")

saveRDS(gwas_cage_bovans_load,
        file = "outputs/gwas_cage_bovans_load_coord.Rds")

saveRDS(gwas_all_bovans_load,
        file = "outputs/gwas_all_bovans_load_coord.Rds")


saveRDS(gwas_pen_lsl_load,
        file = "outputs/gwas_pen_lsl_load_coord.Rds")

saveRDS(gwas_cage_lsl_load,
        file = "outputs/gwas_cage_lsl_load_coord.Rds")

saveRDS(gwas_all_lsl_load,
        file = "outputs/gwas_all_lsl_load_coord.Rds")



saveRDS(gwas_pen_bovans_weight,
        file = "outputs/gwas_pen_bovans_weight_coord.Rds")

saveRDS(gwas_cage_bovans_weight,
        file = "outputs/gwas_cage_bovans_weight_coord.Rds")

saveRDS(gwas_all_bovans_weight,
        file = "outputs/gwas_all_bovans_weight_coord.Rds")


saveRDS(gwas_pen_lsl_weight,
        file = "outputs/gwas_pen_lsl_weight_coord.Rds")

saveRDS(gwas_cage_lsl_weight,
        file = "outputs/gwas_cage_lsl_weight_coord.Rds")

saveRDS(gwas_all_lsl_weight,
        file = "outputs/gwas_all_lsl_weight_coord.Rds")



## Set up chromosome lengths and marker positions for Manhattan plots

chr_lengths <- summarise(group_by(gwas_all_bovans_weight, chr), length = unique(max(ps)))
preferred_order <- c(1:28, 30:31, 33, "Z")

chr_lengths <- chr_lengths[match(preferred_order, chr_lengths$chr),]


gwas_pen_bovans_load$global_pos <- flatten_coordinates(gwas_pen_bovans_load$chr,
                                                       gwas_pen_bovans_load$ps,
                                                       chr_lengths)
gwas_cage_bovans_load$global_pos <- flatten_coordinates(gwas_cage_bovans_load$chr,
                                                        gwas_cage_bovans_load$ps,
                                                        chr_lengths)
gwas_all_bovans_load$global_pos <- flatten_coordinates(gwas_all_bovans_load$chr,
                                                       gwas_all_bovans_load$ps,
                                                       chr_lengths)

gwas_pen_lsl_load$global_pos <- flatten_coordinates(gwas_pen_lsl_load$chr,
                                                    gwas_pen_lsl_load$ps,
                                                    chr_lengths)
gwas_cage_lsl_load$global_pos <- flatten_coordinates(gwas_cage_lsl_load$chr,
                                                     gwas_cage_lsl_load$ps,
                                                     chr_lengths)
gwas_all_lsl_load$global_pos <- flatten_coordinates(gwas_all_lsl_load$chr,
                                                    gwas_all_lsl_load$ps,
                                                    chr_lengths)

gwas_pen_bovans_weight$global_pos <- flatten_coordinates(gwas_pen_bovans_weight$chr,
                                                         gwas_pen_bovans_weight$ps,
                                                         chr_lengths)
gwas_cage_bovans_weight$global_pos <- flatten_coordinates(gwas_cage_bovans_weight$chr,
                                                          gwas_cage_bovans_weight$ps,
                                                          chr_lengths)
gwas_all_bovans_weight$global_pos <- flatten_coordinates(gwas_all_bovans_weight$chr,
                                                         gwas_all_bovans_weight$ps,
                                                         chr_lengths)

gwas_pen_lsl_weight$global_pos <- flatten_coordinates(gwas_pen_lsl_weight$chr,
                                                      gwas_pen_lsl_weight$ps,
                                                      chr_lengths)
gwas_cage_lsl_weight$global_pos <- flatten_coordinates(gwas_cage_lsl_weight$chr,
                                                       gwas_cage_lsl_weight$ps,
                                                       chr_lengths)
gwas_all_lsl_weight$global_pos <- flatten_coordinates(gwas_all_lsl_weight$chr,
                                                      gwas_all_lsl_weight$ps,
                                                      chr_lengths)


chr_breaks <- summarise(group_by(gwas_all_bovans_load, chr),
                        global_pos_break = min(global_pos))

chr_breaks <- chr_breaks[match(preferred_order, chr_breaks$chr),]
chr_breaks$chr_masked <- chr_breaks$chr
chr_breaks$chr_masked[11:32] <- ""


formatting <- list(geom_hline(yintercept = -log10(5e-8),
                              colour = "red",
                              linetype = 2),
                   geom_hline(yintercept = -log10(1e-5),
                              colour = "blue",
                              linetype = 2),
                   theme_bw(),
                   theme(panel.grid = element_blank(),
                         legend.position = "none"),
                   scale_colour_manual(values = c("grey", "black")),
                   xlab(""),
                   ylab(""))



## Main Manhattan plot of breaking strenght and body weight scans

plot_manhattan_cage_bovans_load <- plot_manhattan(gwas_cage_bovans_load,
                                                  p = "p") +
    formatting +
    ylim(0, 10) +
    ggtitle("Bone breaking strength (Bovans, CAGE)")

plot_manhattan_pen_bovans_load <- plot_manhattan(gwas_pen_bovans_load,
                                                 p = "p") +
    formatting +
    ylim(0, 10) +
    ggtitle("Bone breaking strength (Bovans, PEN)")

plot_manhattan_all_bovans_load <- plot_manhattan(gwas_all_bovans_load,
                                                 p = "p") +
    formatting +
    ylim(0, 10) +
    ggtitle("Bone breaking strength (Bovans, JOINT)")



plot_manhattan_cage_lsl_load <- plot_manhattan(gwas_cage_lsl_load,
                                                  p = "p") +
    formatting +
    ylim(0, 10) +
    ggtitle("Bone breaking strength (LSL, CAGE)")

plot_manhattan_pen_lsl_load <- plot_manhattan(gwas_pen_lsl_load,
                                                 p = "p") +
    formatting +
    ylim(0, 10) +
    ggtitle("Bone breaking strength (LSL, PEN)")

plot_manhattan_all_lsl_load <- plot_manhattan(gwas_all_lsl_load,
                                                 p = "p") +
    formatting +
    ylim(0, 10) +
    ggtitle("Bone breaking strength (LSL, JOINT)")

plot_manhattan_load_combined <- ggarrange(plot_manhattan_cage_bovans_load,
                                          plot_manhattan_pen_bovans_load,
                                          plot_manhattan_all_bovans_load,
                                          plot_manhattan_cage_lsl_load,
                                          plot_manhattan_pen_lsl_load,
                                          plot_manhattan_all_lsl_load,
                                          left = "Negative logarithm of p-value",
                                          ncol = 2,
                                          byrow = FALSE)


pdf("figures/plot_manhattan_crosses_separate_load.pdf", width = 10)
print(plot_manhattan_load_combined)
dev.off()




plot_manhattan_cage_bovans_weight <- plot_manhattan(gwas_cage_bovans_weight,
                                                    p = "p") +
    formatting +
    ylim(0, 17) +
    ggtitle("Body weight (Bovans, CAGE)")

plot_manhattan_pen_bovans_weight <- plot_manhattan(gwas_pen_bovans_weight,
                                                   p = "p") +
    formatting +
    ylim(0, 17) +
    ggtitle("Body weight (Bovans, PEN)")

plot_manhattan_all_bovans_weight <- plot_manhattan(gwas_all_bovans_weight,
                                                   p = "p") +
    formatting +
    ylim(0, 17) +
    ggtitle("Body weight (Bovans, JOINT)")


plot_manhattan_cage_lsl_weight <- plot_manhattan(gwas_cage_lsl_weight,
                                                    p = "p") +
    formatting +
    ylim(0, 17) +
    ggtitle("Body weight (LSL, CAGE)")

plot_manhattan_pen_lsl_weight <- plot_manhattan(gwas_pen_lsl_weight,
                                                   p = "p") +
    formatting +
    ylim(0, 17) +
    ggtitle("Body weight (LSL, PEN)")

plot_manhattan_all_lsl_weight <- plot_manhattan(gwas_all_lsl_weight,
                                                   p = "p") +
    formatting +
    ylim(0, 17) +
    ggtitle("Body weight (LSL, JOINT)")


plot_manhattan_weight_combined <- ggarrange(plot_manhattan_cage_bovans_weight,
                                            plot_manhattan_pen_bovans_weight,
                                            plot_manhattan_all_bovans_weight,
                                            plot_manhattan_cage_lsl_weight,
                                            plot_manhattan_pen_lsl_weight,
                                            plot_manhattan_all_lsl_weight,
                                            left = "Negative logarithm of p-value",
                                            ncol = 2,
                                            byrow = FALSE)

pdf("figures/plot_manhattan_crosses_separate_weight.pdf", width = 10)
print(plot_manhattan_weight_combined)
dev.off()



plot_qq_cage_bovans_load <- plot_qq(gwas_cage_bovans_load$p) +
    ggtitle("Bone breaking strength (Bovans, CAGE)")

plot_qq_pen_bovans_load <- plot_qq(gwas_pen_bovans_load$p) +
    ggtitle("Bone breaking strength (Bovans, PEN)")

plot_qq_all_bovans_load <- plot_qq(gwas_all_bovans_load$p) +
    ggtitle("Bone breaking strength (Bovans, JOINT)")


plot_qq_cage_lsl_load <- plot_qq(gwas_cage_lsl_load$p) +
    ggtitle("Bone breaking strength (LSl, CAGE)")

plot_qq_pen_lsl_load <- plot_qq(gwas_pen_lsl_load$p) +
    ggtitle("Bone breaking strength (LSL, PEN)")

plot_qq_all_lsl_load <- plot_qq(gwas_all_lsl_load$p) +
    ggtitle("Bone breaking strength (LSL, JOINT)")


plot_qq_load_combined <- ggarrange(plot_qq_cage_bovans_load,
                                   plot_qq_pen_bovans_load,
                                   plot_qq_all_bovans_load,
                                   plot_qq_cage_lsl_load,
                                   plot_qq_pen_lsl_load,
                                   plot_qq_all_lsl_load,
                                   ncol = 2,
                                   byrow = FALSE)

pdf("figures/plot_qq_crosses_separate_load.pdf", width = 10)
print(plot_qq_load_combined)
dev.off()




plot_qq_cage_bovans_weight <- plot_qq(gwas_cage_bovans_weight$p) +
    ggtitle("Body weight (Bovans, CAGE)")

plot_qq_pen_bovans_weight <- plot_qq(gwas_pen_bovans_weight$p) +
    ggtitle("Body weight (Bovans, PEN)")

plot_qq_all_bovans_weight <- plot_qq(gwas_all_bovans_weight$p) +
    ggtitle("Body weight (Bovans, JOINT)")


plot_qq_cage_lsl_weight <- plot_qq(gwas_cage_lsl_weight$p) +
    ggtitle("Body weight (LSL, CAGE)")

plot_qq_pen_lsl_weight <- plot_qq(gwas_pen_lsl_weight$p) +
    ggtitle("Body weight (LSL, PEN)")

plot_qq_all_lsl_weight <- plot_qq(gwas_all_lsl_weight$p) +
    ggtitle("Body weight (LSL, JOINT)")




plot_qq_weight_combined <- ggarrange(plot_qq_cage_bovans_weight,
                                     plot_qq_pen_bovans_weight,
                                     plot_qq_all_bovans_weight,
                                     plot_qq_cage_lsl_weight,
                                     plot_qq_pen_lsl_weight,
                                     plot_qq_all_lsl_weight,
                                     ncol = 2,
                                     byrow = FALSE)


pdf("figures/plot_qq_crosses_separate_weight.pdf", width = 10)
print(plot_qq_weight_combined)
dev.off()




## Comparison of systems


formatting_comparison <- list(theme_bw(),
                              theme(panel.grid = element_blank()),
                              scale_colour_manual(values = c("blue", "red")),
                              xlab("PEN"),
                              ylab("CAGE"))


gwas_crosses_load <- rbind(transform(gwas_pen_bovans_load,
                                     cross = "Bovans",
                                     scan = "PEN"),
                           transform(gwas_pen_lsl_load,
                                     cross = "LSL", 
                                     scan = "PEN"),
                           transform(gwas_cage_bovans_load,
                                     cross = "Bovans",
                                     scan = "CAGE"),
                           transform(gwas_cage_lsl_load,
                                     cross = "LSL", 
                                     scan = "CAGE"))

crosses_load_comparison <- pivot_wider(gwas_crosses_load,
                                       c("marker_id", "cross", "scan"), 
                                       values_from = c("p", "estimates"),
                                       names_from = c("scan"))


plot_crosses_load_comparison_p <- qplot(x = -log10(p_PEN),
                                        y = -log10(p_CAGE),
                                        colour = cross,
                                        data = filter(crosses_load_comparison,
                                                      p_PEN < 1e-3 |
                                                          p_CAGE < 1e-3)) +
    ggtitle("Bone strength \nNegative logarithm of p-value") +
    formatting_comparison +
    xlim(0, 10) +
    ylim(0, 10)

plot_crosses_load_comparison_beta <- qplot(x = estimates_PEN,
                                           y = estimates_CAGE,
                                           colour = cross,
                                           data = filter(crosses_load_comparison,
                                                         p_PEN < 1e-3 |
                                                             p_CAGE < 1e-3)) +
    geom_hline(yintercept = 0, colour = "red", linetype = 2) +
    geom_vline(xintercept = 0, colour = "red", linetype = 2) +
    ggtitle("Bone strength \nEstimated marker effect") +
    formatting_comparison



gwas_crosses_weight <- rbind(transform(gwas_pen_bovans_weight,
                                       cross = "Bovans",
                                       scan = "PEN"),
                             transform(gwas_pen_lsl_weight,
                                       cross = "LSL", 
                                       scan = "PEN"),
                             transform(gwas_cage_bovans_weight,
                                       cross = "Bovans",
                                       scan = "CAGE"),
                             transform(gwas_cage_lsl_weight,
                                       cross = "LSL", 
                                       scan = "CAGE"))

crosses_weight_comparison <- pivot_wider(gwas_crosses_weight,
                                         c("marker_id", "cross", "scan"), 
                                         values_from = c("p", "estimates"),
                                         names_from = c("scan"))


plot_crosses_weight_comparison_p <- qplot(x = -log10(p_PEN),
                                          y = -log10(p_CAGE),
                                          colour = cross,
                                          data = filter(crosses_weight_comparison,
                                                        p_PEN < 1e-3 |
                                                            p_CAGE < 1e-3)) +
    ggtitle("Body weight \nNegative logarithm of p-value") +
    formatting_comparison 




plot_crosses_weight_comparison_beta <- qplot(x = estimates_PEN,
                                           y = estimates_CAGE,
                                           colour = cross,
                                           data = filter(crosses_weight_comparison,
                                                         p_PEN < 1e-3 |
                                                             p_CAGE < 1e-3)) +
    geom_hline(yintercept = 0, colour = "red", linetype = 2) +
    geom_vline(xintercept = 0, colour = "red", linetype = 2) +
    ggtitle("Body weight \nEstimated marker effect") +
    formatting_comparison


plot_crosses_comparisons <- wrap_plots(list(plot_crosses_load_comparison_p,
                                            plot_crosses_load_comparison_beta,
                                            plot_crosses_weight_comparison_p,
                                            plot_crosses_weight_comparison_beta)) +
    plot_layout(guides = "collect") & theme(legend.title = element_blank())



pdf("figures/plot_gwas_comparison_crosses_separate.pdf", width = 10)
print(plot_crosses_comparisons)
dev.off()








## Create tables

supplementary_load_table <- rbind(transform(gwas_pen_bovans_load[, -7],
                                            scan = "PEN",
                                            cross = "Bovans"),
                                  transform(gwas_cage_bovans_load[, -7],
                                            scan = "CAGE",
                                            cross = "Bovans"),
                                  transform(gwas_all_bovans_load[, -7],
                                            scan = "JOINT",
                                            cross = "Bovans"),
                                  transform(gwas_pen_lsl_load[, -7],
                                            scan = "PEN",
                                            cross = "LSL"),
                                  transform(gwas_cage_lsl_load[, -7],
                                            scan = "CAGE",
                                            cross = "LSL"),
                                  transform(gwas_all_lsl_load[, -7],
                                            scan = "JOINT",
                                            cross = "LSL"))


supplementary_weight_table <- rbind(transform(gwas_pen_bovans_weight[, -7],
                                              scan = "PEN",
                                              cross = "Bovans"),
                                    transform(gwas_cage_bovans_weight[, -7],
                                              scan = "CAGE",
                                              cross = "Bovans"),
                                    transform(gwas_all_bovans_weight[, -7],
                                              scan = "JOINT",
                                              cross = "Bovans"),
                                    transform(gwas_pen_lsl_weight[, -7],
                                              scan = "PEN",
                                              cross = "LSL"),
                                    transform(gwas_cage_lsl_weight[, -7],
                                              scan = "CAGE",
                                              cross = "LSL"),
                                    transform(gwas_all_lsl_weight[, -7],
                                              scan = "JOINT",
                                              cross = "LSL"))




supplementary_suggestive_load <- filter(supplementary_load_table,
                                        p < 1e-5)
supplementary_suggestive_load <-
    supplementary_suggestive_load[order(as.numeric(supplementary_suggestive_load$chr),
                                        supplementary_suggestive_load$ps),]


supplementary_suggestive_weight <- filter(supplementary_weight_table,
                                        p < 1e-6)
supplementary_suggestive_weight <-
    supplementary_suggestive_weight[order(as.numeric(supplementary_suggestive_weight$chr),
                                          supplementary_suggestive_weight$ps),]


write_csv(supplementary_load_table,
          "tables/supplementary_data_bone_strength_crosses_separate_summary_statistics.csv")

write_csv(supplementary_load_table, 
          "tables/supplementary_data_weight_crosses_separate_summary_statistics.csv")

write_csv(supplementary_suggestive_load, 
          "tables/supplementary_table_bone_strength_crosses_separate_hits.csv")

write_csv(supplementary_suggestive_weight, 
          "tables/supplementary_table_weight_crosses_separate_hits.csv")





## Regions to zoom in

plot_hit <- function(hit, formatting) {
    qplot(x = ps/1e6, y = -log10(p), data = hit) +
        formatting
}

weight_bovans_chr4 <- filter(gwas_all_bovans_weight,
                             chr == 4 &
                                 ps > 75748329 - 2e6 &
                                 ps < 75850294 + 2e6)

weight_lsl_chr4 <- filter(gwas_all_lsl_weight,
                          chr == 4 &
                              ps > 75748329 - 2e6 &
                              ps < 75850294 + 2e6)

weight_bovans_chr6 <- filter(gwas_all_bovans_weight,
                             chr == 6 &
                                 ps > 11477631 - 2e6 &
                                 ps < 11477631 + 2e6)

weight_lsl_chr6 <- filter(gwas_all_lsl_weight,
                          chr == 6 &
                              ps > 11477631 - 2e6 &
                              ps < 11477631 + 2e6)

weight_bovans_chr27 <- filter(gwas_all_bovans_weight,
                              chr == 27 &
                                  ps > 6070932 - 2e6 &
                                  ps < 6087051 + 2e6)


weight_lsl_chr27 <- filter(gwas_all_lsl_weight,
                              chr == 27 &
                                  ps > 6070932 - 2e6 &
                                  ps < 6087051 + 2e6)


formatting_weight_hits <- list(geom_hline(yintercept = -log10(5e-8),
                                          colour = "blue",
                                          linetype = 2),
                               theme_bw(),
                               theme(panel.grid = element_blank(),
                                     strip.background = element_blank()),
                               ylim(0, 17),
                               xlab("Position (Mbp)"),
                               ylab("-log(p)"))


plot_weight_hits <- plot_hit(weight_bovans_chr4, formatting_weight_hits) +
    ggtitle("Body weight, chr 4 (Bovans)") +
    plot_hit(weight_lsl_chr4, formatting_weight_hits) +
    ggtitle("Body weight, chr 4 (LSL)") +
    plot_hit(weight_bovans_chr6, formatting_weight_hits) +
    ggtitle("Body weight, chr 6 (Bovans)") + 
    plot_hit(weight_lsl_chr6, formatting_weight_hits) +
    ggtitle("Body weight, chr 6 (LSL)") + 
    plot_hit(weight_bovans_chr27, formatting_weight_hits) +
    ggtitle("Body weight, chr 27 (Bovans)") + 
    plot_hit(weight_lsl_chr27, formatting_weight_hits) +
    ggtitle("Body weight, chr 27 (Bovans)") + 
    plot_layout(ncol = 2)


pdf("figures/plot_weight_hits_crosses_separate.pdf",
    width = 10)
print(plot_weight_hits)
dev.off()




load_chr2_cage <- filter(gwas_cage_bovans_load,
                         chr == 2 & 
                             ps > 125198709 - 2e6 &
                             ps < 125198709 + 2e6)

load_chr2_pen <- filter(gwas_pen_bovans_load,
                        chr == 2 & 
                            ps > 125198709 - 2e6 &
                            ps < 125198709 + 2e6)


load_chr3_cage <- filter(gwas_cage_lsl_load,
                         chr == 3 & 
                             ps > 2430687 - 2e6 &
                             ps < 2430687 + 2e6)

load_chr3_pen <- filter(gwas_pen_lsl_load,
                        chr == 3 & 
                            ps > 2430687 - 2e6 &
                            ps < 2430687 + 2e6)

load_chr10_cage <- filter(gwas_cage_bovans_load,
                          chr == 10 & 
                              ps > 7029216 - 2e6 &
                              ps < 7029216 + 2e6)

load_chr10_pen <- filter(gwas_pen_bovans_load,
                         chr == 10 & 
                             ps > 7029216 - 2e6 &
                             ps < 7029216 + 2e6)


load_chr11_cage <- filter(gwas_cage_lsl_load,
                          chr == 11 & 
                              ps > 18432763 - 2e6 &
                              ps < 18432763 + 2e6)

load_chr11_pen <- filter(gwas_pen_lsl_load,
                         chr == 11 & 
                             ps > 18432763 - 2e6 &
                             ps < 18432763 + 2e6)




formatting_load_hits <- list(geom_hline(yintercept = -log10(1e-5),
                                        colour = "blue",
                                        linetype = 2),
                             theme_bw(),
                             theme(panel.grid = element_blank(),
                                   strip.background = element_blank()),
                             ylim(0, 8),
                             xlab("Position (Mbp)"),
                             ylab("-log(p)"))



plot_load_hits <- plot_hit(load_chr2_cage, formatting_load_hits) +
    ggtitle("Bone strength, chr 2 (Bovans, CAGE)") +
    plot_hit(load_chr2_pen, formatting_load_hits) +
    ggtitle("Bone strength, chr 2 (Bovans, PEN)") +
    plot_hit(load_chr3_cage, formatting_load_hits) +
    ggtitle("Bone strength, chr 3 (LSL, CAGE)") +
    plot_hit(load_chr3_pen, formatting_load_hits) +
    ggtitle("Bone strength, chr 3 (LSL, PEN)") +
    plot_hit(load_chr10_cage, formatting_load_hits) +
    ggtitle("Bone strength, chr 10 (Bovans, CAGE)") +
    plot_hit(load_chr10_pen, formatting_load_hits) +
    ggtitle("Bone strength, chr 10 (Bovans, PEN)") +
    plot_hit(load_chr11_cage, formatting_load_hits) +
    ggtitle("Bone strength, chr 11 (LSL, CAGE)") +
    plot_hit(load_chr11_pen, formatting_load_hits) +
    ggtitle("Bone strength, chr 11 (LSL, PEN)") +
    plot_layout(ncol = 2)

pdf("figures/plot_load_hits_crosses_separate.pdf",
    width = 10)
print(plot_load_hits)
dev.off()





## Conditional GWAS

gwas_weight_conditional <- readRDS("gwas/hglm_gwas_all_weight_conditional.Rds")

gwas_weight_conditional$marker_id <- clean_marker_id(gwas_weight_conditional$marker_id)

gwas_weight_conditional <- inner_join(map, gwas_weight_conditional)

gwas_weight_conditional$global_pos <- flatten_coordinates(gwas_weight_conditional$chr,
                                                          gwas_weight_conditional$ps,
                                                          chr_lengths)




combined_conditional <- rbind(transform(filter(gwas_all_lsl_weight, chr == 4),
                                        scan = "Default"),
                              transform(filter(gwas_weight_conditional, chr == 4),
                                        scan = "Conditional"))

plot_conditional <- qplot(x = ps/1e6, y = -log10(p),
                          colour = scan,
                          data = combined_conditional) +
    formatting +
    ggtitle("Conditional GWAS of body weight on chromosome 4")



pdf("figures/plot_chr4_conditional_lsl.pdf")
print(plot_conditional)
dev.off()
