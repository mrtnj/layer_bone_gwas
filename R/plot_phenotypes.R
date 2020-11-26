
## Exploratory plots of phenotypes

library(broom)
library(dplyr)
library(egg)
library(ggplot2)
library(tidyr)
library(patchwork)
library(readr)


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


## Feed

table(pheno$feed)

table(pheno$feed, pheno$cage.pen)

table(pheno$feed, pheno$cage.pen, pheno$breed)

model_feed <- lm(load_N ~ feed + cage.pen + breed, data = pheno)

print(summary(model_feed))
