
## Exploratory plots of phenotypes

library(broom)
library(dplyr)
library(egg)
library(ggplot2)
library(lme4)
library(tidyr)
library(patchwork)


pheno <- readRDS("outputs/pheno.Rds")

pheno <- filter(pheno, !is.na(cage.pen))

pheno$comb_weight <- pheno$comb_g/pheno$weight
pheno$load_weight <- pheno$load_N/pheno$weight


## Histograms

long <- pivot_longer(pheno,
                     c("load_N", "weight", "comb_g"),
                     names_to = "variable")

plot_histograms <- qplot(x = value, fill = cage.pen, data = long) +
    facet_wrap(variable ~ breed, scale = "free", ncol = 2)


## Histograms of relative phenotypes

long_relative <- pivot_longer(pheno,
                              c("comb_weight", "load_weight"),
                              names_to = "variable")

plot_relative_histograms <- qplot(x = value, fill = cage.pen, data = long_relative) +
    facet_wrap(variable ~ breed, scale = "free", ncol = 2)


## Scatterplots with weight

long_scatter <- pivot_longer(pheno,
                             c("load_N", "comb_g"),
                             names_to = "variable")


plot_scatter <- qplot(x = weight, y = value, colour = cage.pen, data = long_scatter) +
    facet_wrap(variable ~ breed, scale = "free") +
    geom_smooth(method = lm)


## Dotplots of variables with models

model_breed_cagepen <- function(long) {

    long_split <- split(long, long$variable)

    models <- lapply(long_split,
                     function(x) lm(value ~ breed * cage.pen, data = x))
    
    levels <- expand.grid(cage.pen = c("CAGE", "PEN"),
                          breed = c("LSL", "Bovans"),
                          stringsAsFactors = FALSE)
    
    ## Get fits for cells
    
    fit <- lapply(models,
                  predict,
                  newdata = levels,
                  interval = "confidence")
    
    fit_levels <- lapply(fit, cbind, levels)
    
    fits <- Reduce(rbind, fit_levels)
    fits$variable <- rep(names(fit), each = nrow(levels))
    
    
    ## Get comparisons
    
    coef <- lapply(models,
                   tidy,
                   conf.int = TRUE)
    
    coefs <- Reduce(rbind, coef)
    coefs$variable <- rep(names(coef), each = nrow(coef[[1]]))

    list(models = models,
         fits = fits,
         coefs = coefs)
}

models_load_weight <- model_breed_cagepen(long)

fits <- models_load_weight$fits
coefs <- models_load_weight$coefs


plot_dots <- ggplot() +
    geom_jitter(aes(x = cage.pen, colour = breed, y = value),
                alpha = I(0.2),
                data = long) +
    geom_pointrange(aes(x = cage.pen, colour = breed,
                        y = fit, ymin = lwr, ymax = upr),
                    data = fits) +
    facet_wrap(~ variable, scale = "free")




## Regression with body weight

load_split <- split(pheno, list(pheno$breed, pheno$cage.pen))

models_bw <- lapply(load_split,
                    function(split) lm(load_N ~ weight, data = split))

coef_bw <- lapply(models_bw, tidy, conf.int = TRUE)
coefs_bw <- Reduce(rbind, coef_bw)

coefs_bw$model <- rep(names(coef_bw), each = 2)
coefs_bw$breed <- unlist(lapply(strsplit(coefs_bw$model, split = "\\."), "[", 1))
coefs_bw$cage.pen <- unlist(lapply(strsplit(coefs_bw$model, split = "\\."), "[", 2))


## Nicer plots

long_load <- filter(long, variable != "comb_g")
fits_load <- filter(fits, variable != "comb_g")
coefs_load <- filter(coefs,
                     variable != "comb_g" &
                     term != "(Intercept)")

long_load$Variable <- ifelse(long_load$variable == "weight",
                             "Body weight (kg)",
                             "Tibial breaking strength (N)")
fits_load$Variable <- ifelse(fits_load$variable == "weight",
                             "Body weight (kg)",
                             "Tibial breaking strength (N)")
coefs_load$Variable <- ifelse(coefs_load$variable == "weight",
                              "Body weight (kg)",
                              "Tibial breaking strength (N)")

coefs_load$Term <- ""
coefs_load$Term[coefs_load$term == "breedLSL"] <- "Bovans vs LSL"
coefs_load$Term[coefs_load$term == "cage.penPEN"] <- "CAGE vs PEN"
coefs_load$Term[coefs_load$term == "breedLSL:cage.penPEN"] <- "Interaction"

plot_dots_nocomb <- ggplot() +
    geom_jitter(aes(x = cage.pen, colour = breed, y = value),
                alpha = I(0.1),
                data = long_load) +
    geom_pointrange(aes(x = cage.pen, colour = breed,
                        y = fit, ymin = lwr, ymax = upr),
                    position = position_dodge(0.5),
                    data = fits_load) +
    facet_wrap(~ Variable, scale = "free") +
    scale_colour_manual(values = c("blue", "red"),
                        name = "") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.background = element_blank()) +
    xlab("") +
    ylab("")

plot_ests_nocomb <- ggplot() +
    geom_pointrange(aes(x = Term,
                        y = estimate,
                        ymin = conf.low,
                        ymax = conf.high),
                    data = coefs_load) +
    geom_hline(yintercept = 0,
               colour = "red",
               linetype = 2) +
    facet_grid(~Variable, scale = "free") +
    scale_x_discrete(limits = c("Interaction", "Bovans vs LSL", "CAGE vs PEN")) +
    coord_flip() +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.background = element_blank()) +
    xlab("") +
    ylab("") +
    ggtitle("Differences")

plot_trait_combined <- ggarrange(plot_dots_nocomb, plot_ests_nocomb, heights = c(0.7, 0.3))

pdf("figures/plot_housing_breed_difference.pdf")
print(plot_trait_combined)
dev.off()





plot_scatter_nocomb <- qplot(x = weight, y = load_N, colour = breed,
                             alpha = I(0.25),
                             data = pheno) +
    geom_smooth(method = lm, se = FALSE) +
    facet_wrap(~ cage.pen) +
    scale_colour_manual(values = c("blue", "red"),
                        name = "") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.background = element_blank()) +
    xlab("Body weight (kg)") +
    ylab("Tibial breaking strength (N)")


plot_scatter_ests <- ggplot() +
    geom_pointrange(aes(x = cage.pen,
                        colour = breed,
                        y = estimate,
                        ymin = conf.low,
                        ymax = conf.high),
                    position = position_dodge(0.5),
                    data = filter(coefs_bw, term != "(Intercept)")) +
    scale_colour_manual(values = c("blue", "red"),
                        name = "") +
    coord_flip() +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.background = element_blank()) +
    xlab("") +
    ylab("") +
    ggtitle("Regression coefficient (N/kg)")
              


plot_scatter_combined <- ggarrange(plot_scatter_nocomb, plot_scatter_ests, heights = c(0.7, 0.3))


pdf("figures/plot_bw_regression.pdf")
print(plot_scatter_combined)
dev.off()



## pQCT phenotypes

ct_trait_names <- c("TOT_CNT", "TOT_DEN",
                    "TRAB_CNT", "TRAB_DEN",
                    "CRT_CNT", "CRT_DEN",
                    "CRT_THK_C", "OBJECTLEN")

pretty_ct_trait_names <- read_csv("pretty_trait_names_pqct.csv")

ct_ix <- which(colnames(pheno) %in%
               c(paste("ct_mid_", ct_trait_names, sep = ""),
                 paste("ct_distal_", ct_trait_names, sep = "")))

ct_cor <- cor(pheno[, ct_ix], use = "p")


ct_animals <- na.exclude(pheno[, c(1, ct_ix)])$animal_id

ct_pheno <- filter(pheno, animal_id %in% ct_animals)

ct_pca <- prcomp(ct_pheno[, ct_ix], scale = TRUE)

ct_pca_data <- cbind(ct_pheno, ct_pca$x)

long_ct_pca <- pivot_longer(ct_pca_data, PC1:PC9)

ct_pca_loadings <- data.frame(trait = rownames(ct_pca$rotation),
                              ct_pca$rotation[,1:3],
                              stringsAsFactors = FALSE)

ct_pca_loadings_long <- pivot_longer(ct_pca_loadings, -trait,
                                     names_to = "PC")

plot_ct_loadings <- ggplot() +
    geom_bar(aes(x = trait, y = value),
             data = ct_pca_loadings_long,
             stat = "identity") +
    scale_x_discrete(limits = pretty_ct_trait_names$name,
                     labels = pretty_ct_trait_names$pretty_name) +
    facet_wrap(~ PC, ncol = 1) +
    coord_flip() +
    theme_bw() +
    ylab("") +
    xlab("") +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) +
    ggtitle("Principial component loadings")


plot_ct_pcs <- qplot(x = name,
                     y = value,
                     colour = cage.pen,
                     data = long_ct_pca,
                     geom = "boxplot")

plot_ct_load <- qplot(x = load_N,
                      y = value,
                      colour = cage.pen,
                      data = long_ct_pca) +
    facet_wrap(~ name) +
    geom_smooth(method = lm)


## pQCT heatmap

ct_cor_long <- pivot_longer(data.frame(trait1 = rownames(ct_cor),
                                       ct_cor,
                                       stringsAsFactors = FALSE),
                            -trait1,
                            values_to = "correlation",
                            names_to = "trait2")

plot_ct_heatmap <- qplot(x = trait1,
                         y = trait2,
                         fill = correlation, 
                         geom = "tile",
                         data = ct_cor_long) +
    scale_fill_gradient2(limits = c(-1, 1)) +
    scale_x_discrete(limits = pretty_ct_trait_names$name,
                     labels = pretty_ct_trait_names$pretty_name) +
    scale_y_discrete(limits = pretty_ct_trait_names$name,
                     labels = pretty_ct_trait_names$pretty_name) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = -90),
          panel.grid = element_blank()) +
    xlab("") + ylab("") + ggtitle("Pearson correlation between pQCT phenotypes")


ct_pca_variance <- data.frame(var_explained = summary(ct_pca)$importance["Proportion of Variance",])
ct_pca_variance$PC <- sub(rownames(ct_pca_variance),
                          pattern = "PC",
                          replacement = "")


plot_ct_pc_variance <- ggplot() +
    geom_bar(aes(x = PC,
                 y = var_explained),
             data = ct_pca_variance,
             stat = "identity") +
    scale_x_discrete(limits = ct_pca_variance$PC) +
    theme_bw() +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank()) +
    ylim(0, 1) +
    xlab("Principal component") +
    ylab("Proportion variance explained") +
    ggtitle("Variance explained by pQCT principal components")



plot_ct_heatmap_pc_combined <- ggarrange(plot_ct_heatmap, plot_ct_pc_variance,
                                         heights = c(2, 1))

pdf("figures/plot_ct_pca.pdf")
print(plot_ct_heatmap_pc_combined)
dev.off()


pdf("figures/plot_ct_pca_loadings.pdf")
print(plot_ct_loadings)
dev.off()

## FTIR phenotypes

ft_ix <- grep("^ftir_", colnames(pheno))

ft_cor <- cor(pheno[, ft_ix], use = "p")


ft_animals <- na.exclude(pheno[, c(1, ft_ix)])$animal_id

ft_pheno <- filter(pheno, animal_id %in% ft_animals)

ft_pca <- prcomp(ft_pheno[, ft_ix], scale = TRUE)

ft_pca_data <- cbind(ft_pheno, ft_pca$x)

long_ft_pca <- pivot_longer(ft_pca_data, PC1:PC9)

plot_ft_pcs <- qplot(x = name,
                     y = value,
                     colour = cage.pen,
                     data = long_ft_pca,
                     geom = "boxplot")

plot_ft_load <- qplot(x = load_N,
                      y = value,
                      colour = cage.pen,
                      data = long_ft_pca) +
    facet_wrap(~ name,
               scale = "free") +
    geom_smooth(method = lm)



## TGA phenotypes

tga_ix <- grep("MB$|CB$", colnames(pheno))

tga_cor <- cor(pheno[, tga_ix], use = "p")

tga_animals <- na.exclude(pheno[, c(1, tga_ix)])$animal_id

tga_pheno <- filter(pheno, animal_id %in% tga_animals)

tga_pca <- prcomp(tga_pheno[, tga_ix], scale = TRUE)

tga_pca_data <- cbind(tga_pheno, tga_pca$x)

long_tga_pca <- pivot_longer(tga_pca_data, PC1:PC9)

plot_tga_pcs <- qplot(x = name,
                      y = value,
                      colour = cage.pen,
                      data = long_tga_pca,
                      geom = "boxplot")

plot_tga_load <- qplot(x = load_N,
                       y = value,
                       colour = cage.penn,
                       data = long_tga_pca) +
    facet_wrap(~ name,
               scale = "free") +
    geom_smooth(method = lm)


pretty_tga_trait_names <- read_csv("pretty_trait_names_tga.csv")

tga_cor_long <- pivot_longer(data.frame(trait1 = rownames(tga_cor),
                                       tga_cor,
                                       stringsAsFactors = FALSE),
                            -trait1,
                            values_to = "correlation",
                            names_to = "trait2")

plot_tga_heatmap <- qplot(x = trait1,
                         y = trait2,
                         fill = correlation, 
                         geom = "tile",
                         data = tga_cor_long) +
    scale_fill_gradient2(limits = c(-1, 1)) +
    scale_x_discrete(limits = pretty_tga_trait_names$name,
                     labels = pretty_tga_trait_names$pretty_name) +
    scale_y_discrete(limits = pretty_tga_trait_names$name,
                     labels = pretty_tga_trait_names$pretty_name) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = -90),
          panel.grid = element_blank()) +
    xlab("") + ylab("") + ggtitle("Pearson correlation between TGA phenotypes")


tga_pca_loadings <- data.frame(trait = rownames(tga_pca$rotation),
                               tga_pca$rotation[,1:6],
                               stringsAsFactors = FALSE)

tga_pca_loadings_long <- pivot_longer(tga_pca_loadings, -trait,
                                     names_to = "PC")

plot_tga_loadings <- ggplot() +
    geom_bar(aes(x = trait, y = value),
             data = tga_pca_loadings_long,
             stat = "identity") +
    scale_x_discrete(limits = pretty_tga_trait_names$name,
                     labels = pretty_tga_trait_names$pretty_name) +
    facet_wrap(~ PC, ncol = 1) +
    coord_flip() +
    theme_bw() +
    ylab("") +
    xlab("") +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) +
    ggtitle("Principial component loadings")



tga_pca_variance <- data.frame(var_explained = summary(tga_pca)$importance["Proportion of Variance",])
tga_pca_variance$PC <- sub(rownames(tga_pca_variance),
                           pattern = "PC",
                           replacement = "")


plot_tga_pc_variance <- ggplot() +
    geom_bar(aes(x = PC,
                 y = var_explained),
             data = tga_pca_variance,
             stat = "identity") +
    scale_x_discrete(limits = tga_pca_variance$PC) +
    theme_bw() +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank()) +
    ylim(0, 1) +
    xlab("Principal component") +
    ylab("Proportion variance explained") +
    ggtitle("Variance explained by TGA principal components")


plot_tga_heatmap_pc_combined <- ggarrange(plot_tga_heatmap, plot_tga_pc_variance,
                                          heights = c(2, 1))

pdf("figures/plot_tga_pca.pdf")
print(plot_tga_heatmap_pc_combined)
dev.off()


## Figure of differences between systems and association with load


tga_pheno_long <- pivot_longer(tga_pheno[, c(1, 3, 5, tga_ix)],
                               -c(animal_id, cage.pen, breed),
                               names_to = "variable")

models_tga <- model_breed_cagepen(tga_pheno_long)

fits_tga <- models_tga$fits
coefs_tga <- models_tga$coefs


plot_tga_pheno <- ggplot() +
    geom_jitter(aes(colour = breed,
                    y = value,
                    x = cage.pen),
                alpha = I(0.2),
                data = tga_pheno_long) +
    geom_pointrange(aes(colour = breed,
                        y = fit,
                        ymin = lwr,
                        ymax = upr,
                        x = cage.pen),
                    position = position_dodge(0.5),
                    data = fits_tga) +
    facet_wrap(~ variable, scale = "free_y") +
    scale_colour_manual(values = c("blue", "red"),
                        name = "") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.background = element_blank()) +
    xlab("") +
    ylab("")
