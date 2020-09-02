
library(dplyr)
library(egg)
library(ggplot2)
library(patchwork)
library(readr)
library(stringr)
library(tidyr)
library(patchwork)

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


saveRDS(gwas_pen_load,
        file = "outputs/gwas_pen_load_coord.Rds")

saveRDS(gwas_cage_load,
        file = "outputs/gwas_cage_load_coord.Rds")

saveRDS(gwas_all_load,
        file = "outputs/gwas_all_load_coord.Rds")


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


pdf("figures/plot_chr4_conditional_hglm.pdf")
print(plot_conditional)
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




## Candidate regions

suggestive_pen_load <- filter(gwas_pen_load, p < 1e-4)
suggestive_cage_load <- filter(gwas_cage_load, p < 1e-4)
suggestive_all_load <- filter(gwas_all_load, p < 1e-4)

significant_pen_weight <- filter(gwas_pen_weight, p < 5e-8)
significant_cage_weight <- filter(gwas_cage_weight, p < 5e-8)
significant_all_weight <- filter(gwas_all_weight, p < 5e-8)



load_chr4 <- filter(gwas_pen_load,
                    chr == 4 &
                        ps > 60941465 - 2e6 &
                        ps < 60941465 + 2e6)

load_chrZ <- filter(gwas_pen_load,
                    chr == "Z" &
                        ps > 33957226 - 2e6 &
                        ps < 33960416 + 2e6)

load_chr2_locus1 <- filter(gwas_cage_load,
                           chr == 2 &
                               ps > 109944191 - 2e6 &
                               ps < 110288659 + 2e6)

load_chr2_locus2 <- filter(gwas_cage_load,
                           chr == 2 &
                               ps > 118267538 - 2e6 &
                               ps < 125198709 + 2e6)


load_chr3 <- filter(gwas_cage_load,
                    chr == 3 &
                        ps > 91588936 - 2e6 &
                        ps < 91588936 + 2e6)

load_chr11_locus1 <- filter(gwas_cage_load,
                            chr == 11 &
                                ps > 5913130 - 2e6 &
                                ps < 5913130 + 2e6)

load_chr11_locus2 <- filter(gwas_cage_load,
                            chr == 11 &
                                ps > 18476664 - 2e6 &
                                ps < 18476664 + 2e6)




formatting_load_hits <- list(geom_hline(yintercept = -log10(1e-4),
                                        colour = "blue",
                                        linetype = 2),
                             theme_bw(),
                             theme(panel.grid = element_blank(),
                                   strip.background = element_blank()),
                             ylim(0, 8),
                             xlab("Position (Mbp)"),
                             ylab("-log(p)"))

plot_hit <- function(hit, formatting) {
    qplot(x = ps/1e6, y = -log10(p), data = hit) +
        formatting
}

plot_load_hits <- plot_hit(load_chr4, formatting_load_hits) + ggtitle("Bone strength, chr 4 (PEN)") +
    plot_hit(load_chrZ, formatting_load_hits) + ggtitle("Bone strength, chr Z (PEN)") +
    plot_hit(load_chr2_locus1, formatting_load_hits) + ggtitle("Bone strength, chr 2 (CAGE)") +
    plot_hit(load_chr2_locus2, formatting_load_hits) + ggtitle("Bone strength, chr 2 (CAGE)") +
    plot_hit(load_chr3, formatting_load_hits) + ggtitle("Bone strength, chr 3 (CAGE)") +
    plot_hit(load_chr11_locus1, formatting_load_hits) + ggtitle("Bone strength, chr 11 (CAGE)") +
    plot_hit(load_chr11_locus2, formatting_load_hits) + ggtitle("Bone strength, chr 11 (CAGE)") +
    plot_layout(ncol = 2)
                   

pdf("figures/plot_load_hits_hglm.pdf")
print(plot_load_hits)
dev.off()




weight_chr4 <- filter(gwas_all_weight,
                      chr == 4 &
                          ps > 73994772 - 2e6 &
                          ps < 78415388 + 2e6)

weight_chr6 <- filter(gwas_all_weight,
                      chr == 6 &
                          ps > 11403561 - 2e6 &
                          ps < 11558459 + 2e6)

weight_chr27 <- filter(gwas_all_weight,
                      chr == 27 &
                          ps > 6070932 - 2e6 &
                          ps < 6147189 + 2e6)

formatting_weight_hits <- list(geom_hline(yintercept = -log10(5e-8),
                                          colour = "blue",
                                          linetype = 2),
                               theme_bw(),
                               theme(panel.grid = element_blank(),
                                     strip.background = element_blank()),
                               ylim(0, 15),
                               xlab("Position (Mbp)"),
                               ylab("-log(p)"))


plot_weight_hits <- plot_hit(weight_chr4, formatting_weight_hits) + ggtitle("Body weight, chr 4 (JOINT)") +
    plot_hit(weight_chr6, formatting_weight_hits) + ggtitle("Body weight, chr 6 (JOINT)") + 
    plot_hit(weight_chr27, formatting_weight_hits) + ggtitle("Body weight, chr 27 (JOINT)") + 
    plot_layout(ncol = 1)


pdf("figures/plot_weight_hits_hglm.pdf")
print(plot_weight_hits)
dev.off()




## Create tables

supplementary_load_table <- rbind(transform(gwas_pen_load[, -7], scan = "PEN"),
                                  transform(gwas_cage_load[, -7], scan = "CAGE"),
                                  transform(gwas_all_load[, -7], scan = "JOINT"))

supplementary_weight_table <- rbind(transform(gwas_pen_weight[, -7], scan = "PEN"),
                                    transform(gwas_cage_weight[, -7], scan = "CAGE"),
                                    transform(gwas_all_weight[, -7], scan = "JOINT"))


supplementary_suggestive_load <- filter(supplementary_load_table,
                                        p < 1e-4)
supplementary_suggestive_load <-
    supplementary_suggestive_load[order(as.numeric(supplementary_suggestive_load$chr),
                                                                     supplementary_suggestive_load$ps),]


supplementary_suggestive_weight <- filter(supplementary_weight_table,
                                        p < 5e-8)
supplementary_suggestive_weight <-
    supplementary_suggestive_weight[order(as.numeric(supplementary_suggestive_weight$chr),
                                          supplementary_suggestive_weight$ps),]


write_csv(supplementary_load_table,
          "tables/supplementary_data_bone_strength_summary_statistics.csv")

write_csv(supplementary_load_table, 
          "tables/supplementary_data_weight_summary_statistics.csv")

write_csv(supplementary_suggestive_load, 
          "tables/supplementary_table_bone_strength_hits.csv")

write_csv(supplementary_suggestive_weight, 
          "tables/supplementary_table_weight_hits.csv")




## GALLO analysis of regions

library(GALLO)


qtl_annotation <- import_gff_gtf("annotation/chickenqtldb2020-06_16.txt", "gff")

qtl_bed <- data.frame(chr = paste("chr", qtl_annotation$chr, sep = ""),
                      start = qtl_annotation$start_pos - 1,
                      end = qtl_annotation$end_pos,
                      qtlid = 1:nrow(qtl_annotation),
                      stringsAsFactors = FALSE)

write.table(qtl_bed,
            file = "temp_qtl_bed.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

## Lift with UCSC liftOver

qtl_bed_lifted <- read_tsv("annotation/chickenqtldb_liftover.bed",
                           col_names = FALSE)
colnames(qtl_bed_lifted) <- c("chr_lifted", "start_lifted", "end_lifted", "qtlid")

qtl_annotation$qtlid <- 1:nrow(qtl_annotation)

qtl_annotation_lifted <- inner_join(qtl_bed_lifted,
                                    qtl_annotation)

qtl_annotation_lifted$start_pos <- qtl_annotation_lifted$start_lifted
qtl_annotation_lifted$end_pos <- qtl_annotation_lifted$end_lifted

qtl_annotation_lifted <- qtl_annotation_lifted[5:10]

weight_gallo <- supplementary_suggestive_weight[, c("chr", "ps")]
colnames(weight_gallo) <- c("CHR", "BP")


load_gallo <- supplementary_suggestive_load[, c("chr", "ps")]
colnames(load_gallo) <- c("CHR", "BP")


gallo_overlap_weight <- find_genes_qtls_around_markers(db_file = as.data.frame(qtl_annotation),
                                                       marker_file = as.data.frame(weight_gallo),
                                                       marker = "snp",
                                                       method = "qtl")

gallo_overlap_load <- find_genes_qtls_around_markers(db_file = as.data.frame(qtl_annotation),
                                                     marker_file = as.data.frame(load_gallo),
                                                     marker = "snp",
                                                     method = "qtl")



gallo_enrichment_weight <- qtl_enrich(qtl_db = as.data.frame(qtl_annotation_lifted),
                                      qtl_file = as.data.frame(qtl_overlap_weight),
                                      qtl_type = "Name",
                                      enrich_type = "genome",
                                      chr.subset = NULL,
                                      padj = "holm")

gallo_enrichment_load <- qtl_enrich(qtl_db = as.data.frame(qtl_annotation_lifted),
                                    qtl_file = as.data.frame(gallo_overlap_load),
                                    qtl_type = "Name",
                                    enrich_type = "genome",
                                    chr.subset = NULL,
                                    padj = "holm")
           




GALLO::QTLenrich_plot(gallo_enrichment_weight, x = "QTL", pval = "adj.pval")

GALLO::QTLenrich_plot(gallo_enrichment_load, x = "QTL", pval = "adj.pval")




