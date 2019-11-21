
## Exploratory plots of phenotypes

library(dplyr)
library(ggplot2)
library(lme4)
library(tidyr)


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

long_split <- split(long, long$variable)

models <- lapply(long_split,
                 function(x) lm(value ~ breed * cage.pen, data = x))

levels <- expand.grid(cage.pen = c("CAGE", "PEN"),
                      breed = c("LSL", "Bovans"),
                      stringsAsFactors = FALSE)

fit <- lapply(models,
              predict,
              newdata = levels,
              interval = "confidence")

fit_levels <- lapply(fit, cbind, levels)

fits <- Reduce(rbind, fit_levels)
fits$variable <- rep(names(fit), each = nrow(levels))

plot_dots <- ggplot() +
    geom_jitter(aes(x = cage.pen, colour = breed, y = value),
                alpha = I(0.2),
                data = long) +
    geom_pointrange(aes(x = cage.pen, colour = breed,
                        y = fit, ymin = lwr, ymax = upr),
                    data = fits) +
    facet_wrap(~ variable, scale = "free")



## Model with group effect

model <- lmer(load_N ~ breed * cage.pen + (1 | group),
              data = pheno)
