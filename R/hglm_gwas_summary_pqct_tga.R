
## Summarise results of bone phenotype scans of traits with significant h2

library(dplyr)
library(egg)
library(ggplot2)
library(patchwork)
library(readr)
library(stringr)
library(tidyr)
library(patchwork)

source("R/gwas_helper_functions.R")


gwas_pen <- readRDS("gwas/hglm_gwas_pen_bone_phenotypes.Rds")
gwas_cage <- readRDS("gwas/hglm_gwas_cage_bone_phenotypes.Rds")
gwas_all <- readRDS("gwas/hglm_gwas_all_bone_phenotypes.Rds")



clean_marker_id <- function(marker_id) {
    sub(marker_id, pattern = "_[ATGC]", replacement = "")
}

gwas_pen$marker_id <- clean_marker_id(gwas_pen$marker_id)
gwas_cage$marker_id <- clean_marker_id(gwas_cage$marker_id)
gwas_all$marker_id <- clean_marker_id(gwas_all$marker_id) 



map <- read_delim("gwas/all.map",
                  delim = " ",
                  col_types = "ccnn",
                  col_names = FALSE)[,-3]

colnames(map) <- c("chr", "marker_id", "ps")

gwas_pen <- inner_join(map, gwas_pen)
gwas_cage <- inner_join(map, gwas_cage)
gwas_all <- inner_join(map, gwas_all)


pretty_trait_names <- rbind(data.frame(name = c("ct_pc1", "ct_pc2", "ct_pc3"),
                                       pretty_name = c("pQCT PC1 'high density, thickness, content'",
                                                             "pQCT PC2 'long bone length'",
                                                             "pQCT PC3 'low cortical density'")),
                            read_csv("pretty_trait_names_tga.csv"))
colnames(pretty_trait_names)[1] <- "trait"


gwas_pen <- inner_join(pretty_trait_names, gwas_pen)
gwas_cage <- inner_join(pretty_trait_names, gwas_cage)
gwas_all <- inner_join(pretty_trait_names, gwas_all)


saveRDS(gwas_pen,
        file = "outputs/gwas_pen_bone_phenotypes_coord.Rds")

saveRDS(gwas_cage,
        file = "outputs/gwas_cage_bone_phenotypes_coord.Rds")

saveRDS(gwas_all,
        file = "outputs/gwas_all_bone_phenotypes_coord.Rds")


## Set up chromosome lengths and marker positions for Manhattan plots

chr_lengths <- summarise(group_by(gwas_all, chr), length = unique(max(ps)))
preferred_order <- c(1:28, 30:31, 33, "Z")

chr_lengths <- chr_lengths[match(preferred_order, chr_lengths$chr),]

gwas_pen$global_pos <- flatten_coordinates(gwas_pen$chr, gwas_pen$ps, chr_lengths)
gwas_cage$global_pos <- flatten_coordinates(gwas_cage$chr, gwas_cage$ps, chr_lengths)
gwas_all$global_pos <- flatten_coordinates(gwas_all$chr, gwas_all$ps, chr_lengths)

chr_breaks <- summarise(group_by(gwas_all, chr),
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


pen_traits <- unique(gwas_pen$trait)
cage_traits <- unique(gwas_cage$trait)
all_traits <- unique(gwas_all$trait)

pen_only_traits <- setdiff(pen_traits, cage_traits)
cage_only_traits <- setdiff(cage_traits, pen_traits)



## One scan only

manhattans_pen_only <- lapply(pen_only_traits,
                              function(trait) {
                                  plot_manhattan(gwas_pen[gwas_pen$trait == trait,],
                                                 p = "p") + 
                                      formatting +
                                      ylim(0, 10) +
                                      ggtitle(paste(pretty_trait_names$pretty_name[pretty_trait_names$trait == trait],
                                                    "(PEN)"))
                                  
                              })

plot_manhattan_pen_only <- do.call(ggarrange, manhattans_pen_only)



manhattans_cage_only <- lapply(cage_only_traits,
                              function(trait) {
                                  plot_manhattan(gwas_cage[gwas_cage$trait == trait,],
                                                 p = "p") +
                                      formatting +
                                      ylim(0, 10) +
                                      ggtitle(paste(pretty_trait_names$pretty_name[pretty_trait_names$trait == trait],
                                                    "(CAGE)"))
                              })

plot_manhattan_cage_only <- do.call(ggarrange, manhattans_cage_only)

plot_manhattan_one_only <- do.call(ggarrange, c(manhattans_cage_only, manhattans_pen_only, ncol = 2))

pdf("figures/plot_manhattan_pqct_tga_one_only.pdf", width = 10)
print(plot_manhattan_one_only)
dev.off()


## Multi-scan traits


manhattans_cage_multi <- lapply(all_traits,
                                function(trait) {
                                    plot_manhattan(gwas_cage[gwas_cage$trait == trait,],
                                                   p = "p") +
                                        formatting +
                                        ggtitle(paste(pretty_trait_names$pretty_name[pretty_trait_names$trait == trait],
                                                      "(CAGE)"))
                                })

manhattans_pen_multi <- lapply(all_traits,
                                function(trait) {
                                    plot_manhattan(gwas_pen[gwas_pen$trait == trait,],
                                                   p = "p") +
                                        formatting +
                                        ggtitle(paste(pretty_trait_names$pretty_name[pretty_trait_names$trait == trait],
                                                      "(PEN)"))
                                })

manhattans_all_multi <- lapply(all_traits,
                                function(trait) {
                                    plot_manhattan(gwas_all[gwas_all$trait == trait,],
                                                   p = "p") +
                                        formatting +
                                        ggtitle(paste(pretty_trait_names$pretty_name[pretty_trait_names$trait == trait],
                                                      "(JOINT)"))
                                })


plot_manhattan_multi_pqct <- ggarrange(manhattans_cage_multi[[1]], manhattans_pen_multi[[1]], manhattans_all_multi[[1]],
                                       manhattans_cage_multi[[2]], manhattans_pen_multi[[2]], manhattans_all_multi[[2]],
                                       byrow = FALSE)

pdf("figures/plot_manhattan_pqct_all.pdf", width = 10)
print(plot_manhattan_multi_pqct)
dev.off()

plot_manhattan_multi_tga <- ggarrange(manhattans_cage_multi[[3]] + ylim(0, 10),
                                      manhattans_pen_multi[[3]] + ylim(0, 10),
                                      manhattans_all_multi[[3]] + ylim(0, 10),
                                      manhattans_cage_multi[[4]] + ylim(0, 10),
                                      manhattans_pen_multi[[4]] + ylim(0, 10),
                                      manhattans_all_multi[[4]] + ylim(0, 10),
                                      byrow = FALSE)

pdf("figures/plot_manhattan_tga_all.pdf", width = 10)
print(plot_manhattan_multi_tga)
dev.off()

