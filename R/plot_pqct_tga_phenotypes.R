
## Models and plots of content, density and composition phenotypes

library(ggplot2)
library(multcomp)
library(tidyr)
library(patchwork)
library(readr)

pheno <- readRDS("outputs/pheno.Rds")

pheno <- filter(pheno, !is.na(cage.pen))



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


tga_pheno_long <- pivot_longer(tga_pheno[, c(1, 3, 5, 7, tga_ix)],
                               -c(animal_id, cage.pen, breed, weight),
                               names_to = "variable")


model_breed_cagepen <- function(long) {
    
    long_split <- split(long, long$variable)
    
    models <- lapply(long_split,
                     function(x) {
                         lm(value ~ breed * cage.pen, data = x)
                     })
    
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
    
    
    ## Get contrasts
    cage_minus_pen_within_bovan <- matrix(c(0, 0, -1, 0), 1)
    cage_minus_pen_within_lsl <- matrix(c(0, 0, -1, -1), 1)
    
    contr <- lapply(models, function(model) {
        bovan_contrast <- tidy(confint(glht(model, cage_minus_pen_within_bovan)))
        lsl_contrast <- tidy(confint(glht(model, cage_minus_pen_within_lsl)))
        
        rbind(transform(bovan_contrast, contrast = "within Bovans"),
              transform(lsl_contrast, contrast = "within LSL"))
        })
    
    contr_df <- Reduce(rbind, contr)
    contr_df$variable <- rep(names(contr), each = nrow(contr[[1]]))
    
    list(models = models,
         fits = fits,
         coefs = coefs,
         contr = contr_df)
}




models_tga <- model_breed_cagepen_bw(tga_pheno_long)

fits_tga <- models_tga$fits
coefs_tga <- models_tga$coefs
contr_tga <- models_tga$contr

coefs_tga$Term <- ""
coefs_tga$Term[coefs_load$term == "breedLSL"] <- "Bovans vs LSL"
coefs_tga$Term[coefs_load$term == "cage.penPEN"] <- "CAGE vs PEN"
coefs_tga$Term[coefs_load$term == "breedLSL:cage.penPEN"] <- "Interaction"

coefs_tga <- inner_join(coefs_tga, pretty_tga_trait_names,
                        by = c("variable" = "name"))
coefs_tga$pretty_name <- factor(coefs_tga$pretty_name,
                                levels = pretty_tga_trait_names$pretty_name[order(pretty_tga_trait_names$name)])

fits_tga <- inner_join(fits_tga, pretty_tga_trait_names,
                       by = c("variable" = "name"))
fits_tga$pretty_name <- factor(fits_tga$pretty_name,
                               levels = pretty_tga_trait_names$pretty_name[order(pretty_tga_trait_names$name)])

tga_pheno_long_pretty <- inner_join(tga_pheno_long, pretty_tga_trait_names,
                                    by = c("variable" = "name"))
tga_pheno_long_pretty$pretty_name <- factor(tga_pheno_long_pretty$pretty_name,
                                            levels = pretty_tga_trait_names$pretty_name[order(pretty_tga_trait_names$name)])

plot_tga_pheno <- ggplot() +
    geom_jitter(aes(colour = breed,
                    y = value,
                    x = cage.pen),
                alpha = I(0.2),
                data = tga_pheno_long_pretty) +
    geom_pointrange(aes(colour = breed,
                        y = fit,
                        ymin = lwr,
                        ymax = upr,
                        x = cage.pen),
                    position = position_dodge(0.5),
                    data = fits_tga) +
    facet_wrap(~ pretty_name, scale = "free_y", ncol = 3) +
    scale_colour_manual(values = c("blue", "red"),
                        name = "") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.background = element_blank()) +
    xlab("") +
    ylab("")


plot_tga_coef <- qplot(x = Term,
                       y = estimate,
                       ymin = conf.low,
                       ymax = conf.high,
                       data = filter(coefs_tga, term != "(Intercept)"), geom = "pointrange") +
    facet_wrap(~ pretty_name, scale = "free", ncol = 2) +
    geom_hline(yintercept = 0,
               colour = "red",
               linetype = 2) +
    scale_x_discrete(limits = c("Interaction", "Bovans vs LSL", "CAGE vs PEN")) +
    coord_flip() +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.background = element_blank()) +
    xlab("") +
    ylab("") 


plot_tga_contr <- qplot(x = contrast,
                        y = estimate,
                        ymin = conf.low,
                        ymax = conf.high,
                        geom = "pointrange",
                        data = contr_tga) +
    facet_wrap(~ variable, scale = "free", ncol = 2) +
    geom_hline(yintercept = 0,
               colour = "red",
               linetype = 2) +
    coord_flip() +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.background = element_blank()) +
    xlab("") +
    ylab("") 

pdf("figures/plot_tga_pheno.pdf")
print(plot_tga_pheno)
dev.off()

pdf("figures/plot_tga_coefs.pdf")
print(plot_tga_coef)
dev.off()


