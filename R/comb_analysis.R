
## Analyse comb data

library(broom)
library(dplyr)
library(egg)
library(genetics)
library(ggplot2)
library(readr)
library(tidyr)


source("R/gwas_helper_functions.R")



pheno <- readRDS("outputs/pheno.Rds")

pheno <- filter(pheno, !is.na(cage.pen) & !is.na(comb_g))

pheno$comb_weight <- pheno$comb_g/pheno$weight
pheno$comb_residual <- residuals(lm(comb_g ~ weight, pheno, na.action = na.omit))


## Histograms

long <- pivot_longer(pheno,
                     c("weight", "comb_g", "comb_weight", "comb_residual"),
                     names_to = "variable")

plot_histograms <- qplot(x = value, fill = cage.pen, data = long) +
    facet_wrap(variable ~ breed, scale = "free", ncol = 2)


## Scatterplots with weight

plot_scatter <- qplot(x = weight, y = comb_g, colour = cage.pen, data = pheno) +
    facet_wrap( ~ breed, scale = "free") +
    geom_smooth(method = lm) +
    scale_colour_manual(values = c("blue", "red"), name = "") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.background = element_blank()) +
    ggtitle("Relationship between comb mass and body mass") +
    ylab("Comb mass (g)") +
    xlab("Body mass (kg)")

## Regression with body weight

comb_split <- split(pheno, list(pheno$breed, pheno$cage.pen))

models_bw <- lapply(comb_split,
                    function(split) lm(comb_g ~ weight, data = split))

coef_bw <- lapply(models_bw, tidy, conf.int = TRUE)
coefs_bw <- Reduce(rbind, coef_bw)

coefs_bw$model <- rep(names(coef_bw), each = 2)
coefs_bw$breed <- unlist(lapply(strsplit(coefs_bw$model, split = "\\."), "[", 1))
coefs_bw$cage.pen <- unlist(lapply(strsplit(coefs_bw$model, split = "\\."), "[", 2))


plot_slope <- ggplot() +
    geom_pointrange(aes(x = cage.pen,
                        colour = breed,
                        y = estimate,
                        ymin = conf.low,
                        ymax = conf.high),
                    position = position_dodge(0.5),
                    data = filter(coefs_bw, term != "(Intercept)")) +
    scale_colour_manual(values = c("blue", "red"), name = "") +
    scale_y_continuous(breaks = seq(from = 0, to = 12, by = 2)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.background = element_blank()) +
    coord_flip() +
    ylab("Comb mass (g) / Body mass (kg)") +
    xlab("Housing system")

model_slope_interaction <- lm(comb_g ~ weight * cage.pen, pheno)

plot_scatter_slope <- ggarrange(plot_scatter,
                                plot_slope,
                                heights = c(0.75, 0.25),
                                labels = c("A", "B"))


pdf("figures/comb_scatterplot.pdf")
print(plot_scatter_slope)
dev.off()



## Dotplots of variables with models

model_residual <- lm(comb_residual ~ breed * cage.pen, data = pheno)

levels <- expand.grid(cage.pen = c("CAGE", "PEN"),
                      breed = c("LSL", "Bovans"),
                      stringsAsFactors = FALSE)

## Get fits for cells

fit <- predict(model_residual,
               newdata = levels,
               interval = "confidence")

fit_levels <- cbind(fit, levels)


## Get comparisons

coef_residual <- tidy(model_residual,
                      conf.int = TRUE)

coef_residual$Term <- ""
coef_residual$Term[coef_residual$term == "breedLSL"] <- "LSL vs Bovans"
coef_residual$Term[coef_residual$term == "cage.penPEN"] <- "PEN vs CAGE"
coef_residual$Term[coef_residual$term == "breedLSL:cage.penPEN"] <- "Interaction"


plot_dots <- ggplot() +
    geom_jitter(aes(x = cage.pen, colour = breed, y = comb_residual),
                alpha = I(0.2),
                data = pheno) +
    geom_pointrange(aes(x = cage.pen, colour = breed,
                        y = fit, ymin = lwr, ymax = upr),
                    data = fit_levels,
                    position = position_dodge(0.5)) +
    scale_colour_manual(values = c("blue", "red"), name = "") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.background = element_blank()) +
    ggtitle("Differences in residual comb mass bewteen housing systems\nand crossbreds") +
    ylab("Residual comb mass (g)") +
    xlab("Housing system")



plot_ests <- ggplot() +
    geom_pointrange(aes(x = Term,
                        y = estimate,
                        ymin = conf.low,
                        ymax = conf.high),
                    position = position_dodge(0.5),
                    data = filter(coef_residual, term != "(Intercept)")) +
    scale_colour_manual(values = c("blue", "red"),
                        name = "") +
    scale_y_continuous(breaks = seq(from = -6, to = 3, by = 1)) +
    scale_x_discrete(limits = c("Interaction", "LSL vs Bovans", "PEN vs CAGE")) +
    geom_hline(yintercept = 0, colour = "red", linetype = 2) +
    coord_flip() +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.background = element_blank()) +
    xlab("") +
    ylab("") +
    ggtitle("Estimated differences (g)")


plot_dots_ests <- ggarrange(plot_dots,
                            plot_ests,
                            heights = c(0.75, 0.25),
                            labels = c("A", "B"))


pdf("figures/comb_dotplot.pdf")
print(plot_dots_ests)
dev.off()



## Investigating the min, max, and range of groups

group_range <- function(pheno) {
    summarise(group_by(pheno, group, cage.pen, breed),
              mean = mean(comb_residual),
              min = min(comb_residual),
              max = max(comb_residual),
              range = max - min)
}

group_size <- as.data.frame(table(pheno$group, pheno$breed),
                            stringsAsFactors = FALSE)
colnames(group_size) <- c("group", "breed", "size")
group_size <- filter(group_size, size != 0)
group_size$cage.pen <- ifelse(group_size$group > 2000,
                              "CAGE",
                              "PEN")



sim_groups <- function() {
    do(group_by(group_size, group, cage.pen, breed), {
        random_inds <- sample(1:nrow(pheno), .$size, replace = TRUE)
        data.frame(comb_residual = pheno$comb_residual[random_inds])
    })
}

sim_group_size <- replicate(100, sim_groups(), simplify = FALSE)



real_group_range <- group_range(pheno)
sim_group_range <- lapply(sim_group_size, group_range)


long_real_range <- pivot_longer(real_group_range, -c("group", "cage.pen", "breed"))
long_sim_range <- pivot_longer(sim_group_range[[1]], -c("group", "cage.pen", "breed"))


long_real_range$Name <- ""
long_real_range$Name[long_real_range$name == "min"] <- "Group minimum"
long_real_range$Name[long_real_range$name == "mean"] <- "Group mean"
long_real_range$Name[long_real_range$name == "max"] <- "Group maximum"
long_real_range$Name <- factor(long_real_range$Name,
                               levels = c("Group minimum", "Group mean", "Group maximum"))


plot_real_range <- qplot(x = value, fill = cage.pen, data = long_real_range) +
    facet_wrap(breed ~ name, scale = "free")

plot_sim_range <- qplot(x = value, fill = cage.pen, data = long_sim_range) +
    facet_wrap(breed ~ name, scale = "free")
    

plot_sim_vs_real <- ggarrange(plot_real_range, plot_sim_range)


sim_max <- unlist(lapply(sim_group_range,
                         function(x) mean(x$max)))


plot_group_ranges_bars <- qplot(x = factor(group),
                                y = mean,
                                ymin = min,
                                ymax = max,
                                data = real_group_range,
                                colour = cage.pen,
                                geom = "pointrange") +
    scale_colour_manual(values = c("grey", "black"), name = "") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    xlab("Group") +
    ylab("Range and mean of residual comb mass") +
    facet_wrap(~breed)


plot_max_min_histograms <- qplot(x = value, fill = cage.pen,
                                 data = filter(long_real_range,
                                               name %in% c("max", "mean", "min"))) +
    facet_wrap(breed ~ Name) +
    scale_fill_manual(values = c("grey", "black"), name = "") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.background = element_blank()) +
    xlab("Residual comb mass (g)") +
    ylab("Number of groups in bin")


plot_ranges_combined <- ggarrange(plot_group_ranges_bars,
                                  plot_max_min_histograms,
                                  heights = c(0.6, 0.4))


pdf("figures/comb_group_ranges.pdf")
print(plot_ranges_combined)
dev.off()



## GWAS

gwas_pen <- read_tsv("gwas/pen_comb_g/output/pen_comb_g.assoc.txt",
                     col_types = "ccnnccnnnnnnnnn")
gwas_cage <- read_tsv("gwas/cage_comb_g/output/cage_comb_g.assoc.txt",
                      col_types = "ccnnccnnnnnnnnn")
gwas_all <- read_tsv("gwas/all_comb_g/output/all_comb_g.assoc.txt",
                     col_types = "ccnnccnnnnnnnnn")



chr_lengths <- summarise(group_by(gwas_pen, chr), length = unique(max(ps)))
preferred_order <- c(1:28, 30:31, 33, "Un_NW_020110160v1", "Un_NW_020110165v1", "Z")

chr_lengths <- chr_lengths[match(preferred_order, chr_lengths$chr),]


gwas_pen$global_pos <- flatten_coordinates(gwas_pen$chr, gwas_pen$ps, chr_lengths)
gwas_cage$global_pos <- flatten_coordinates(gwas_cage$chr, gwas_cage$ps, chr_lengths)
gwas_all$global_pos <- flatten_coordinates(gwas_all$chr, gwas_all$ps, chr_lengths)

chr_breaks <- summarise(group_by(gwas_pen, chr),
                        global_pos_break = min(global_pos))

chr_breaks <- chr_breaks[match(preferred_order, chr_breaks$chr),]
chr_breaks$chr_masked <- chr_breaks$chr
chr_breaks$chr_masked[11:33] <- ""


formatting_manhattan <- list(theme_bw(),
                             theme(panel.grid = element_blank(),
                                   legend.position = "none"),
                             scale_colour_manual(values = c("black", "grey")),
                             geom_hline(yintercept = 5, colour = "blue", linetype = 2),
                             ylab(""),
                             xlab(""),
                             ylim(0, 8))

plot_manhattan_pen <- plot_manhattan(gwas_pen) +
    formatting_manhattan +
    ggtitle("Comb mass (PEN)")
plot_manhattan_cage <- plot_manhattan(gwas_cage) +
    formatting_manhattan +
    ggtitle("Comb mass (CAGE)")
plot_manhattan_all <- plot_manhattan(gwas_all) +
    formatting_manhattan +
    ggtitle("Comb mass (JOINT)")

plot_manhattan_combined <- ggarrange(plot_manhattan_cage,
                                     plot_manhattan_pen,
                                     plot_manhattan_all,
                                     left = "Negative logarithm of p-value",
                                     bottom = "Chromosome")

pdf("figures/comb_manhattan.pdf")
print(plot_manhattan_combined)
dev.off()


plot_qq_pen <- plot_qq(gwas_pen$p_wald) +
    ggtitle("Comb mass (PEN)")
plot_qq_cage <- plot_qq(gwas_cage$p_wald) +
    ggtitle("Comb mass (CAGE)")
plot_qq_all <- plot_qq(gwas_all$p_wald) +
    ggtitle("Comb mass (JOINT)")

plot_qq_combined <- ggarrange(plot_qq_cage,
                              plot_qq_pen,
                              plot_qq_all)

pdf("figures/comb_qq.pdf")
print(plot_qq_combined)
dev.off()


suggestive_hit <- filter(gwas_pen, p_wald < 1e-5)

focal_marker <- suggestive_hit$rs
focal_chr <- suggestive_hit$chr

region_pen <- filter(gwas_pen,
                     chr == focal_chr &
                     ps > suggestive_hit$ps - 2e6 &
                     ps < suggestive_hit$ps + 2e6)

region_cage <- filter(gwas_cage,
                      chr == focal_chr &
                      ps > suggestive_hit$ps - 2e6 &
                      ps < suggestive_hit$ps + 2e6)

region_all <- filter(gwas_all,
                     chr == focal_chr &
                     ps > suggestive_hit$ps - 2e6 &
                     ps < suggestive_hit$ps + 2e6)


## LD

geno <- readRDS("outputs/geno.Rds")

get_ld <- function(geno, region, focal_marker) {

    geno <- geno[,colnames(geno) %in% region$rs]
    
    focal_geno <- as.genotype(as.data.frame(geno)[, focal_marker])
    
    r2 <- numeric(ncol(geno))
    
    for(col_ix in 1:ncol(geno)) {
        r2[col_ix] <- LD(focal_geno,
                         as.genotype(as.data.frame(geno)[, col_ix]))$"R^2"
    }
    names(r2) <- colnames(geno)

    r2
}

region_pen$ld <- get_ld(geno[geno$individual >= 1000,],
                        region_pen,
                        focal_marker)

region_cage$ld <- get_ld(geno[geno$individual < 1000,],
                        region_cage,
                        focal_marker)

region_all$ld <- get_ld(geno,
                        region_all,
                        focal_marker)



## Combined plot

formatting_hits <- list(theme_bw(),
                        theme(panel.grid = element_blank(),
                              legend.position = "none"),
                        scale_colour_gradient(low = "blue", high = "red", name = "LD"),
                        geom_hline(yintercept = 5, colour = "blue", linetype = 2),
                        ylab(""),
                        xlab(""),
                        ylim(0, 8))

plot_hit_pen <- qplot(x = ps/1e6, y = -log10(p_wald), colour = ld,
                      data = region_pen) +
    formatting_hits + theme(legend.position = "right") +
    ggtitle("Chromosome 15 locus (PEN)")

plot_hit_cage <- qplot(x = ps/1e6, y = -log10(p_wald), colour = ld,
                      data = region_cage) +
    formatting_hits + 
    ggtitle("Chromosome 15 locus (CAGE)")

plot_hit_all <- qplot(x = ps/1e6, y = -log10(p_wald), colour = ld,
                      data = region_all) +
    formatting_hits +
    ggtitle("Chromosome 15 locus (JOINT)")

plot_hit_combined <- ggarrange(plot_hit_cage,
                               plot_hit_pen,
                               plot_hit_all,
                               left = "Negative logarithm of p-value",
                               bottom = "Position (Mbp)")



pdf("figures/comb_chr15_locus.pdf")
print(plot_hit_combined)
dev.off()
