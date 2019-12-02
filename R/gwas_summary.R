
## Read GWAS results and make graphs

library(dplyr)
library(egg)
library(ggplot2)
library(qqman)
library(readr)
library(stringr)


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


saveRDS(gwas,
        file = "outputs/gwas.Rds")



## Plots

chr_lengths <- summarise(group_by(gwas, chr), length = unique(max(ps)))
preferred_order <- c(1:28, 30:31, 33, "Un_NW_020110160v1", "Un_NW_020110165v1", "Z")

chr_lengths <- chr_lengths[match(preferred_order, chr_lengths$chr),]

flatten_coordinates <- function(chr,
                                pos,
                                chr_lengths) {
    pos_flat <- pos
    offset <- 0
 
    for (chr_ix in 1:nrow(chr_lengths)) {
        on_chr <- chr == chr_lengths$chr[chr_ix]
        pos_flat[on_chr] <- pos[on_chr] + offset
        offset <- offset + chr_lengths$length[chr_ix]
    }
 
    pos_flat
}


gwas$global_pos <- flatten_coordinates(gwas$chr, gwas$ps, chr_lengths)


chr_breaks <- summarise(group_by(gwas, chr),
                        global_pos_break = min(global_pos))

chr_breaks <- chr_breaks[match(preferred_order, chr_breaks$chr),]
chr_breaks$chr_masked <- chr_breaks$chr
chr_breaks$chr_masked[11:33] <- ""



## Manhattan plots


plot_manhattan <- function(data) {
    chr_numbers <- as.numeric(factor(data$chr,
                                     levels = unique(data$chr)))
    data$chr_indicator <- factor(chr_numbers %% 2)

    qplot(x = global_pos,
          y = -log10(p_wald),
          colour = chr_indicator,
          size = I(0.5),
          data = data) +
        scale_x_continuous(breaks = chr_breaks$global_pos_break,
                           labels = chr_breaks$chr_masked)
}


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
                   ylab(""),
                   ylim(0, 10))

plot_manhattan_cage_load <- plot_manhattan(filter(gwas, scan_name == "cage_load_adj")) +
    formatting +
    ggtitle("Bone breaking strength (CAGE)")
plot_manhattan_pen_load <- plot_manhattan(filter(gwas, scan_name == "pen_load_adj")) +
    formatting +
        ggtitle("Bone breaking strength (PEN)")
plot_manhattan_all_load <- plot_manhattan(filter(gwas, scan_name == "all_load_adj")) +
    formatting +
    ggtitle("Bone breaking strength (JOINT)")


plot_manhattan_cage_weight <- plot_manhattan(filter(gwas, scan_name == "cage_weight")) +
    formatting +
    ggtitle("Body weight (CAGE)")
plot_manhattan_pen_weight <- plot_manhattan(filter(gwas, scan_name == "pen_weight")) +
    formatting +
    ggtitle("Body weight (PEN)")
plot_manhattan_all_weight <- plot_manhattan(filter(gwas, scan_name == "all_weight")) +
    formatting +
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

plot_qq <- function(p) {

    Observed <- -log10(sort(p, decreasing = FALSE))
    Expected <- -log10(ppoints(length(p)))

    qplot(x = Expected,
          y = Observed) +
        geom_abline(intercept = 0, slope = 1) +
        theme_bw() +
        theme(panel.grid = element_blank())
}


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
