
## Read GWAS results and make graphs

library(dplyr)
library(egg)
library(ggplot2)
library(qqman)
library(readr)
library(stringr)
library(tidyr)




## Collate results

files <- system("ls gwas/*/output/*assoc.txt", intern = TRUE)

scan_name <- str_match(files, "/([a-z_]+)\\.assoc\\.txt")[,2]

results <- lapply(files,
                  read_tsv,
                  col_types = "ccnnccnnnnnnnnn")

for (file_ix in 1:length(results)) {
    results[[file_ix]]$scan_name <- scan_name[file_ix]
    results[[file_ix]]$numeric_chr <- as.numeric(results[[file_ix]]$chr)
}

## qq(results[[8]]$p_wald)
## manhattan(na.exclude(results[[5]]), chr = "numeric_chr", p = "p_wald", bp = "ps", snp = "rs")



gwas <- Reduce(rbind, results)

gwas$scan_name_pretty <- gwas$scan_name
gwas$scan_name_pretty[gwas$scan_name_pretty == "cage_weight"] <- "Body weight (CAGE)"
gwas$scan_name_pretty[gwas$scan_name_pretty == "pen_weight"] <- "Body weight (PEN)"
gwas$scan_name_pretty[gwas$scan_name_pretty == "all_weight"] <- "Body weight (JOINT)"
gwas$scan_name_pretty[gwas$scan_name_pretty == "all_load_adj"] <-
    "Tibial breaking strength (JOINT)"
gwas$scan_name_pretty[gwas$scan_name_pretty == "pen_load_adj"] <-
    "Tibial breaking strength (PEN)"

gwas$scan_name_pretty <- factor(gwas$scan_name_pretty,
                                levels=c("Body weight (CAGE)",
                                         "Body weight (PEN)",
                                         "Body weight (JOINT)",
                                         "Tibial breaking strength (PEN)",
                                         "Tibial breaking strength (JOINT)"))

gwas$chr_numeric <- as.numeric(gwas$chr)

saveRDS(gwas,
        file = "outputs/gwas.Rds")



## Plots

chr_lengths <- summarise(group_by(gwas, chr), length = unique(max(ps)))
preferred_order <- c(1:28, 30:31, 33, "Un_NW_020110160v1", "Un_NW_020110165v1", "Z")

chr_lengths <- chr_lengths[match(preferred_order, chr_lengths$chr),]

gwas$global_pos <- flatten_coordinates(gwas$chr, gwas$ps, chr_lengths)


chr_breaks <- summarise(group_by(gwas, chr),
                        global_pos_break = min(global_pos))

chr_breaks <- chr_breaks[match(preferred_order, chr_breaks$chr),]
chr_breaks$chr_masked <- chr_breaks$chr
chr_breaks$chr_masked[11:33] <- ""



## Manhattan plots

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

plot_manhattan_cage_load <- plot_manhattan(filter(gwas, scan_name == "cage_load_adj")) +
    formatting +
    ylim(0, 10) +
    ggtitle("Bone breaking strength (CAGE)")
plot_manhattan_pen_load <- plot_manhattan(filter(gwas, scan_name == "pen_load_adj")) +
    formatting +
    ylim(0, 10) +
    ggtitle("Bone breaking strength (PEN)")
plot_manhattan_all_load <- plot_manhattan(filter(gwas, scan_name == "all_load_adj")) +
    formatting +
    ylim(0, 10) +
    ggtitle("Bone breaking strength (JOINT)")


plot_manhattan_cage_weight <- plot_manhattan(filter(gwas, scan_name == "cage_weight")) +
    formatting +
    ylim(0, 17) +
    ggtitle("Body weight (CAGE)")
plot_manhattan_pen_weight <- plot_manhattan(filter(gwas, scan_name == "pen_weight")) +
    formatting +
    ylim(0, 17) +
    ggtitle("Body weight (PEN)")
plot_manhattan_all_weight <- plot_manhattan(filter(gwas, scan_name == "all_weight")) +
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



pdf("figures/plot_manhattans.pdf")
print(plot_manhattan_combined)
dev.off()



## QQ-plots

plot_qq_cage_load <- plot_qq(filter(gwas, scan_name == "cage_load_adj")$p_wald) +
    ggtitle("Bone breaking strength (CAGE)")

plot_qq_pen_load <- plot_qq(filter(gwas, scan_name == "pen_load_adj")$p_wald) +
    ggtitle("Bone breaking strength (PEN)")

plot_qq_all_load <- plot_qq(filter(gwas, scan_name == "all_load_adj")$p_wald) +
    ggtitle("Bone breaking strength (JOINT)")



plot_qq_cage_weight <- plot_qq(filter(gwas, scan_name == "cage_weight")$p_wald) +
    ggtitle("Body weight (CAGE)")

plot_qq_pen_weight <- plot_qq(filter(gwas, scan_name == "pen_weight")$p_wald) +
    ggtitle("Body weight (PEN)")

plot_qq_all_weight <- plot_qq(filter(gwas, scan_name == "all_weight")$p_wald) +
    ggtitle("Body weight (JOINT)")



plot_qq_combined <- ggarrange(plot_qq_cage_load,
                              plot_qq_pen_load,
                              plot_qq_all_load,
                              plot_qq_cage_weight,
                              plot_qq_pen_weight,
                              plot_qq_all_weight,
                              ncol = 2,
                              byrow = FALSE)



pdf("figures/plot_qq.pdf")
print(plot_qq_combined)
dev.off()


## Candidate regions

significant <- filter(gwas, p_wald < 5e-8)

significant_regions <- summarise(group_by(significant, chr, scan_name),
                                 start = min(ps),
                                 end = max(ps),
                                 width = end - start,
                                 mean_af = mean(af),
                                 mean_beta = mean(beta),
                                 mean_p = mean(p_wald),
                                 n_markers = n())



## Weight hits

formatting_weight_hits <- list(geom_hline(yintercept = -log10(5e-8),
                                          colour = "red",
                                          linetype = 2),
                               theme_bw(),
                               theme(panel.grid = element_blank(),
                                     strip.background = element_blank()),
                               ylim(0, 17),
                               xlab(""),
                               ylab(""))

chr4_weight <- filter(gwas, chr == 4 &
                            ps >= 73994772 - 1e6 &
                            ps <= 75850294 + 1e6 &
                            scan_name %in% c("all_weight", "pen_weight", "cage_weight"))

plot_chr4 <- qplot(x = ps/1e6, y = -log10(p_wald), data = chr4_weight) +
    facet_wrap(~ scan_name_pretty) +
    ggtitle("Chromosome 4 locus") +
    formatting_weight_hits

chr6_weight <- filter(gwas, chr == 6 &
                            ps >= 11403561 - 1e6 &
                            ps <= 11588664 + 1e6 &
                            scan_name %in% c("all_weight", "pen_weight", "cage_weight"))

plot_chr6 <- qplot(x = ps/1e6, y = -log10(p_wald), data = chr6_weight) +
    facet_wrap(~ scan_name_pretty) +
    ggtitle("Chromosome 6 locus") +
    formatting_weight_hits

chr27_weight <- filter(gwas, chr == 27 &
                             ps >= 6070932 - 1e6 &
                             ps <= 6147189 + 1e6 &
                             scan_name %in% c("all_weight", "pen_weight", "cage_weight"))

plot_chr27 <- qplot(x = ps/1e6, y = -log10(p_wald), data = chr27_weight) +
    facet_wrap(~ scan_name_pretty) +
    ggtitle("Chromosome 27 locus") +
    formatting_weight_hits



plot_weight_hits <- ggarrange(plot_chr4,
                              plot_chr6,
                              plot_chr27,
                              left = "Negative logarithm of p-value",
                              bottom = "Position (Mbp)")

pdf("figures/plot_gwas_weight_hits.pdf")
print(plot_weight_hits)
dev.off()


## Load hits

suggestive <- filter(gwas, p_wald < 1e-4)

suggestive_regions <- summarise(group_by(suggestive, chr, scan_name),
                                start = min(ps),
                                end = max(ps),
                                width = end - start,
                                mean_af = mean(af),
                                mean_beta = mean(beta),
                                mean_p = mean(p_wald),
                                n_markers = n())



suggestive_load <- filter(suggestive_regions,
                          scan_name %in% c("pen_load_adj", "cage_load_adj", "all_load_adj"))


suggestive_load_list <- vector(mode = "list",
                               length = nrow(suggestive_load))

for (region_ix in 1:nrow(suggestive_load)) {
    this_region <- gwas[gwas$scan_name == suggestive_load$scan_name[region_ix] &
                        gwas$chr == suggestive_load$chr[region_ix] &
                        gwas$ps >= suggestive_load$start[region_ix] - 1e6 &
                        gwas$ps <= suggestive_load$start[region_ix] + 1e6,]
    suggestive_load_list[[region_ix]] <- this_region
}


load_candidates <- Reduce(rbind, suggestive_load_list)

formatting_load_hits <- list(geom_hline(yintercept = -log10(1e-4),
                                          colour = "blue",
                                          linetype = 2),
                               theme_bw(),
                               theme(panel.grid = element_blank(),
                                     strip.background = element_blank()),
                               ylim(0, 10),
                               xlab("Position (Mbp)"),
                               ylab("Negative logarithm of p-value"))

plot_load_hits <- qplot(x = ps/1e6, y = -log10(p_wald), data = load_candidates) +
    facet_wrap(~ chr_numeric + scan_name_pretty, scale = "free_x") +
    formatting_load_hits
    


pdf("figures/plot_gwas_load_suggestive.pdf")
print(plot_load_hits)
dev.off()



## Correlation between scans

gwas_load_systems <- filter(gwas,
                            scan_name %in% c("pen_load_adj", "cage_load_adj"))

load_comparison <- pivot_wider(gwas_load_systems[, c("rs", "p_wald", "beta", "scan_name")],
                               values_from = c("beta", "p_wald"),
                               names_from = "scan_name")


gwas_weight_systems <- filter(gwas,
                            scan_name %in% c("pen_weight", "cage_weight"))

weight_comparison <- pivot_wider(gwas_weight_systems[, c("rs", "p_wald", "beta", "scan_name")],
                                 values_from = c("beta", "p_wald"),
                                 names_from = "scan_name")

formatting_comparison <- list(theme_bw(),
                              theme(panel.grid = element_blank()),
                              xlab("CAGE"),
                              ylab("PEN"))


plot_load_comparison_p <- qplot(x = -log10(p_wald_cage_load_adj),
                                y = -log10(p_wald_pen_load_adj),
                                data = filter(load_comparison,
                                              p_wald_cage_load_adj < 1e-3 |
                                              p_wald_pen_load_adj < 1e-3)) +
    ggtitle("Tibial breaking strength \nNegative logarithm of p-value") +
    formatting_comparison +
    xlim(0, 10) +
    ylim(0, 10)

plot_load_comparison_beta <- qplot(x = beta_cage_load_adj,
                                   y = beta_pen_load_adj,
                                   data = filter(load_comparison,
                                                 p_wald_cage_load_adj < 1e-3 |
                                                 p_wald_pen_load_adj < 1e-3)) +
    geom_hline(yintercept = 0, colour = "red", linetype = 2) +
    geom_vline(xintercept = 0, colour = "red", linetype = 2) +
    ggtitle("Tibial breaking strength \nEstimated marker effect") +
    formatting_comparison


plot_weight_comparison_p <- qplot(x = -log10(p_wald_cage_weight),
                                  y = -log10(p_wald_pen_weight),
                                  data = filter(weight_comparison,
                                                p_wald_cage_weight < 1e-3 |
                                                p_wald_pen_weight < 1e-3)) +
    ggtitle("Body weight \nNegative logarithm of p-value") +
    formatting_comparison +
    xlim(0, 10) +
    ylim(0, 10)

plot_weight_comparison_beta <- qplot(x = beta_cage_weight,
                                     y = beta_pen_weight,
                                     data = filter(weight_comparison,
                                                   p_wald_cage_weight < 1e-3 |
                                                   p_wald_pen_weight < 1e-3)) +
    geom_hline(yintercept = 0, colour = "red", linetype = 2) +
    geom_vline(xintercept = 0, colour = "red", linetype = 2) +
    ggtitle("Body weight \nEstimated marker effect") +
    formatting_comparison


plot_comparisons <- ggarrange(plot_load_comparison_p,
                              plot_load_comparison_beta,
                              plot_weight_comparison_p,
                              plot_weight_comparison_beta)



pdf("figures/plot_gwas_comparison.pdf")
print(plot_comparisons)
dev.off()

cor.test(filter(weight_comparison,
                p_wald_cage_weight < 1e-3 |
                p_wald_pen_weight < 1e-3)$beta_cage_weight,
         filter(weight_comparison,
                p_wald_cage_weight < 1e-3 |
                p_wald_pen_weight < 1e-3)$beta_pen_weight)

cor.test(filter(load_comparison,
                p_wald_cage_load_adj < 1e-3 |
                p_wald_pen_load_adj < 1e-3)$beta_cage_load_adj,
         filter(load_comparison,
                p_wald_cage_load_adj < 1e-3 |
                p_wald_pen_load_adj < 1e-3)$beta_pen_load_adj)


## Tables

significant_table <- significant[, c("scan_name", "rs", "chr", "ps", "p_wald", "beta")]


write.csv(significant_table,
          file = "tables/gwas_significant.csv",
          quote = FALSE,
          row.names = FALSE)


suggestive_table <- suggestive_load[, c("scan_name", "rs", "chr", "ps", "p_wald", "beta")]


write.csv(candidate_table,
          file = "tables/gwas_suggestive.csv",
          quote = FALSE,
          row.names = FALSE)
