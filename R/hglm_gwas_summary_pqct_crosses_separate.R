
## Summarise gwas for pQCT with separate crosses

library(dplyr)
library(egg)
library(ggplot2)
library(patchwork)
library(purrr)
library(readr)
library(stringr)
library(tidyr)
library(patchwork)

source("R/gwas_helper_functions.R")


gwas_pen_bovans <- readRDS("gwas/hglm_gwas_pen_bovans_ct.Rds")
gwas_cage_bovans <- readRDS("gwas/hglm_gwas_cage_bovans_ct.Rds")
gwas_all_bovans <- readRDS("gwas/hglm_gwas_all_bovans_ct.Rds")

gwas_pen_lsl <- readRDS("gwas/hglm_gwas_pen_lsl_ct.Rds")
gwas_cage_lsl <- readRDS("gwas/hglm_gwas_cage_lsl_ct.Rds")
gwas_all_lsl <- readRDS("gwas/hglm_gwas_all_lsl_ct.Rds")



clean_marker_id <- function(marker_id) {
    sub(marker_id, pattern = "_[ATGC]", replacement = "")
}

gwas_pen_bovans$marker_id <- clean_marker_id(gwas_pen_bovans$marker_id)
gwas_cage_bovans$marker_id <- clean_marker_id(gwas_cage_bovans$marker_id)
gwas_all_bovans$marker_id <- clean_marker_id(gwas_all_bovans$marker_id) 

gwas_pen_lsl$marker_id <- clean_marker_id(gwas_pen_lsl$marker_id)
gwas_cage_lsl$marker_id <- clean_marker_id(gwas_cage_lsl$marker_id)
gwas_all_lsl$marker_id <- clean_marker_id(gwas_all_lsl$marker_id) 



map <- read_delim("gwas/all.map",
                  delim = " ",
                  col_types = "ccnn",
                  col_names = FALSE)[,-3]

colnames(map) <- c("chr", "marker_id", "ps")

gwas_pen_bovans <- inner_join(map, gwas_pen_bovans)
gwas_cage_bovans <- inner_join(map, gwas_cage_bovans)
gwas_all_bovans <- inner_join(map, gwas_all_bovans)

gwas_pen_lsl <- inner_join(map, gwas_pen_lsl)
gwas_cage_lsl <- inner_join(map, gwas_cage_lsl)
gwas_all_lsl <- inner_join(map, gwas_all_lsl)



pretty_trait_names <- rbind(data.frame(name = c("ct_pc1", "ct_pc2", "ct_pc3"),
                                       pretty_name = c("QCT PC1 'high density, thickness, content'",
                                                             "QCT PC2 'long bone length'",
                                                             "QCT PC3 'low cortical density'")),
                            read_csv("pretty_trait_names_tga.csv")[,1:2])
colnames(pretty_trait_names)[1] <- "trait"


gwas_pen_bovans <- inner_join(pretty_trait_names, gwas_pen_bovans)
gwas_cage_bovans <- inner_join(pretty_trait_names, gwas_cage_bovans)
gwas_all_bovans <- inner_join(pretty_trait_names, gwas_all_bovans)

gwas_pen_lsl <- inner_join(pretty_trait_names, gwas_pen_lsl)
gwas_cage_lsl <- inner_join(pretty_trait_names, gwas_cage_lsl)
gwas_all_lsl <- inner_join(pretty_trait_names, gwas_all_lsl)



saveRDS(gwas_pen_bovans,
         file = "outputs/gwas_pen_bovans_bone_phenotypes_coord.Rds")

saveRDS(gwas_cage_bovans,
        file = "outputs/gwas_cage_bovans_bone_phenotypes_coord.Rds")

saveRDS(gwas_all_bovans,
        file = "outputs/gwas_all_bovans_bone_phenotypes_coord.Rds")


saveRDS(gwas_pen_lsl,
        file = "outputs/gwas_pen_lsl_bone_phenotypes_coord.Rds")

saveRDS(gwas_cage_lsl,
        file = "outputs/gwas_cage_lsl_bone_phenotypes_coord.Rds")

saveRDS(gwas_all_lsl,
        file = "outputs/gwas_all_lsl_bone_phenotypes_coord.Rds")



## Set up chromosome lengths and marker positions for Manhattan plots

chr_lengths <- summarise(group_by(gwas_all_bovans, chr), length = unique(max(ps)))
preferred_order <- c(1:28, 30:31, 33, "Z")

chr_lengths <- chr_lengths[match(preferred_order, chr_lengths$chr),]

gwas_pen_bovans$global_pos <- flatten_coordinates(gwas_pen_bovans$chr, gwas_pen_bovans$ps, chr_lengths)
gwas_cage_bovans$global_pos <- flatten_coordinates(gwas_cage_bovans$chr, gwas_cage_bovans$ps, chr_lengths)
gwas_all_bovans$global_pos <- flatten_coordinates(gwas_all_bovans$chr, gwas_all_bovans$ps, chr_lengths)

gwas_pen_lsl$global_pos <- flatten_coordinates(gwas_pen_lsl$chr, gwas_pen_lsl$ps, chr_lengths)
gwas_cage_lsl$global_pos <- flatten_coordinates(gwas_cage_lsl$chr, gwas_cage_lsl$ps, chr_lengths)
gwas_all_lsl$global_pos <- flatten_coordinates(gwas_all_lsl$chr, gwas_all_lsl$ps, chr_lengths)



chr_breaks <- summarise(group_by(gwas_all_bovans, chr),
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
                   ylab(""),
                   ylim(0, 15))



gwas_pen_bovans$global_pos <- flatten_coordinates(gwas_pen_bovans$chr, gwas_pen_bovans$ps, chr_lengths)
gwas_cage_bovans$global_pos <- flatten_coordinates(gwas_cage_bovans$chr, gwas_cage_bovans$ps, chr_lengths)
gwas_all_bovans$global_pos <- flatten_coordinates(gwas_all_bovans$chr, gwas_all_bovans$ps, chr_lengths)

gwas_pen_lsl$global_pos <- flatten_coordinates(gwas_pen_lsl$chr, gwas_pen_lsl$ps, chr_lengths)
gwas_cage_lsl$global_pos <- flatten_coordinates(gwas_cage_lsl$chr, gwas_cage_lsl$ps, chr_lengths)
gwas_all_lsl$global_pos <- flatten_coordinates(gwas_all_lsl$chr, gwas_all_lsl$ps, chr_lengths)



make_trait_manhattan <- function(trait, gwas, scan_label) {
    plot_manhattan(gwas[gwas$trait == trait,],
                   p = "p") +
        formatting +
        ggtitle(str_wrap(paste(pretty_trait_names$pretty_name[pretty_trait_names$trait == trait],
                               scan_label),
                         width = 50))
}

manhattans_bovans_pen <- map(c("ct_pc1", "ct_pc2", "ct_pc3"),
                             make_trait_manhattan,
                             gwas = gwas_pen_bovans,
                             scan_label = "(Bovans, PEN)")

manhattans_bovans_cage <- map(c("ct_pc1", "ct_pc2", "ct_pc3"),
                              make_trait_manhattan,
                              gwas = gwas_cage_bovans,
                              scan_label = "(Bovans, CAGE)")

manhattans_bovans_all <- map(c("ct_pc1", "ct_pc2", "ct_pc3"),
                             make_trait_manhattan,
                             gwas = gwas_all_bovans,
                             scan_label = "(Bovans, JOINT)")



manhattans_lsl_pen <- map(c("ct_pc1", "ct_pc2", "ct_pc3"),
                             make_trait_manhattan,
                             gwas = gwas_pen_lsl,
                             scan_label = "(LSL, PEN)")

manhattans_lsl_cage <- map(c("ct_pc1", "ct_pc2", "ct_pc3"),
                              make_trait_manhattan,
                              gwas = gwas_cage_lsl,
                              scan_label = "(LSL, CAGE)")

manhattans_lsl_all <- map(c("ct_pc1", "ct_pc2", "ct_pc3"),
                             make_trait_manhattan,
                             gwas = gwas_all_lsl,
                             scan_label = "(LSL, JOINT)")


plot_manhattan_pc1 <- manhattans_bovans_cage[[1]] + ylim(0, 10) +
    manhattans_lsl_cage[[1]] + ylim(0, 10) + 
    manhattans_bovans_pen[[1]] + ylim(0, 10) +
    manhattans_lsl_pen[[1]] + ylim(0, 10) +
    manhattans_bovans_all[[1]] + ylim(0, 10) + 
    manhattans_lsl_all[[1]] + ylim(0, 10) + 
    plot_layout(ncol = 2)


plot_manhattan_pc2 <- manhattans_bovans_cage[[2]] + ylim(0, 15) +
    manhattans_lsl_cage[[2]] + ylim(0, 15) + 
    manhattans_bovans_pen[[2]] + ylim(0, 15) +
    manhattans_lsl_pen[[2]] + ylim(0, 15) +
    manhattans_bovans_all[[2]] + ylim(0, 15) + 
    manhattans_lsl_all[[2]] + ylim(0, 15) + 
    plot_layout(ncol = 2)


plot_manhattan_pc3 <- manhattans_bovans_cage[[3]] + ylim(0, 10) +
    manhattans_lsl_cage[[3]] + ylim(0, 10) + 
    manhattans_bovans_pen[[3]] + ylim(0, 10) +
    manhattans_lsl_pen[[3]] + ylim(0, 10) +
    manhattans_bovans_all[[3]] + ylim(0, 10) + 
    manhattans_lsl_all[[3]] + ylim(0, 10) + 
    plot_layout(ncol = 2)




pdf("figures/plot_manhattan_pqct1_crosses_separate.pdf", width = 10)
print(plot_manhattan_pc1)
dev.off()


pdf("figures/plot_manhattan_pqct2_crosses_separate.pdf", width = 10)
print(plot_manhattan_pc2)
dev.off()


pdf("figures/plot_manhattan_pqct3_crosses_separate.pdf", width = 10)
print(plot_manhattan_pc3)
dev.off()




## Create tables

supplementary_table <- rbind(transform(gwas_pen_bovans[, -c(1, 9)],
                                       scan = "PEN",
                                       cross = "Bovans"),
                             transform(gwas_cage_bovans[, -c(1, 9)],
                                       scan = "CAGE",
                                       cross = "Bovans"),
                             transform(gwas_all_bovans[, -c(1, 9)],
                                       scan = "JOINT",
                                       cross = "Bovans"),
                             transform(gwas_pen_lsl[, -c(1, 9)],
                                       scan = "PEN",
                                       cross = "LSL"),
                             transform(gwas_cage_lsl[, -c(1, 9)],
                                       scan = "CAGE",
                                       cross = "LSL"),
                             transform(gwas_all_lsl[, -c(1, 9)],
                                       scan = "JOINT",
                                       cross = "LSL"))

colnames(supplementary_table)[1] <- "trait"


supplementary_suggestive <- filter(supplementary_table,
                                   p < 1e-4)
supplementary_suggestive <-
    supplementary_suggestive[order(as.numeric(supplementary_suggestive$chr),
                                   supplementary_suggestive$ps),]


write_csv(supplementary_table,
          "tables/supplementary_data_qct_crosses_separate_summary_statistics.csv")

write_csv(supplementary_suggestive, 
          "tables/supplementary_table_qct_crosses_separate_hits.csv")
