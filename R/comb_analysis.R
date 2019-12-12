
## Analyse comb data

library(broom)
library(dplyr)
library(egg)
library(ggplot2)
library(tidyr)


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
                                heights = c(0.75, 0.25))


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
coef_residual$Term[coef_residual$term == "breedLSL"] <- "Bovans vs LSL"
coef_residual$Term[coef_residual$term == "cage.penPEN"] <- "CAGE vs PEN"
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
    scale_x_discrete(limits = c("Interaction", "Bovans vs LSL", "CAGE vs PEN")) +
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
                            heights = c(0.75, 0.25))


pdf("figures/comb_dotplot.pdf")
print(plot_dots_ests)
dev.off()



##


group_range <- summarise(group_by(pheno, group, cage.pen), mean = mean(comb_residual), min = min(comb_residual), max = max(comb_residual), range = max - min)
