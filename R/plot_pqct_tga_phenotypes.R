
## Models and plots of content, density and composition phenotypes

library(broom)
library(dplyr)
library(egg)
library(ggplot2)
library(multcomp)
library(tidyr)
library(patchwork)
library(purrr)
library(readr)

source("R/trait_modelling_functions.R")

pheno <- readRDS("outputs/pheno.Rds")

pheno <- filter(pheno, !is.na(cage.pen))

write.csv(pheno[!grepl("^ftir_", colnames(pheno))],
          file = "tables/supplementary_data3_phenotypes.csv",
          quote = FALSE,
          row.names = FALSE)



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
    xlab("") + ylab("") + ggtitle("Pearson correlation between QCT phenotypes")


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
    ggtitle("Variance explained by QCT principal components")



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


pretty_tga_trait_names <- read_csv("pretty_trait_names_tga.csv")[,1:2]

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


## Differences between systems: CT PCs


ct_model_data_long <- pivot_longer(ct_pca_data[, c("animal_id",
                                                   "cage.pen",
                                                   "breed",
                                                   "PC1",
                                                   "PC2",
                                                   "PC3")],
                                   -c(animal_id, cage.pen, breed),
                                   names_to = "variable")
                                   


models_ct <- model_breed_cagepen(ct_model_data_long)

contr_ct <- models_ct$contr

contr_ct$pretty_name <- rep(c("QCT PC1 'high density, thickness, content'",
                              "QCT PC2 'long bone length'",
                              "QCT PC3 'low cortical density'"),
                            each = 2)



plot_ct_contr <- qplot(x = contrast,
                       y = estimate,
                       ymin = conf.low,
                       ymax = conf.high,
                       geom = "pointrange",
                       data = contr_ct) +
    facet_wrap(~ variable, scale = "free", ncol = 1) +
    geom_hline(yintercept = 0,
               colour = "red",
               linetype = 2) +
    coord_flip() +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.background = element_blank()) +
    xlab("") +
    ylab("") 


## Figure of differences between systems and association with load

tga_pheno_long <- pivot_longer(tga_pheno[, c(1, 3, 5, tga_ix)],
                               -c(animal_id, cage.pen, breed),
                               names_to = "variable")


models_tga <- model_breed_cagepen(tga_pheno_long)

fits_tga <- models_tga$fits
coefs_tga <- models_tga$coefs
contr_tga <- models_tga$contr
    
coefs_tga$Term <- ""
coefs_tga$Term[coefs_tga$term == "breedLSL"] <- "Bovans vs LSL"
coefs_tga$Term[coefs_tga$term == "cage.penPEN"] <- "CAGE vs PEN"
coefs_tga$Term[coefs_tga$term == "breedLSL:cage.penPEN"] <- "Interaction"

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


tga_contr_pretty <- inner_join(contr_tga, pretty_tga_trait_names,
                               by = c("variable" = "name"))
tga_contr_pretty$pretty_name <- factor(tga_contr_pretty$pretty_name,
                                            levels = pretty_tga_trait_names$pretty_name[order(pretty_tga_trait_names$name)])


fits_ct <- models_ct$fits
fits_ct$pretty_name <- rep(c("QCT PC1 'high density, thickness, content'",
                             "QCT PC2 'long bone length'",
                             "QCT PC3 'low cortical density'"),
                            each = 2)


combined_tga_ct_fits <- rbind(fits_tga,
                              fits_ct)

# combined_tga_ct_contr$pretty_name <- factor(combined_tga_ct_contr$pretty_name,
#                                             levels = c(pretty_tga_trait_names$pretty_name[c(8, 1, 9, 2, 10, 3, 11,
#                                                                                             4, 12, 5, 13, 6, 14, 7)],
#                                                        unique(contr_ct$pretty_name)))


## Estimates for CT and TGA together

plot_tga_estimates <- ggplot() +
    geom_pointrange(aes(colour = breed,
                        y = fit,
                        ymin = lwr,
                        ymax = upr,
                        x = cage.pen),
                    position = position_dodge(0.5),
                    data = combined_tga_ct_fits) +
    facet_wrap(~ pretty_name, scale = "free_y", ncol = 4) +
    scale_colour_manual(values = c("blue", "red"),
                        name = "") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.background = element_blank()) +
    xlab("") +
    ylab("")


## Same plot, reduced to main TGA phenotypes

main_tga_ct_phenotypes <- c("Phosphates_over_OM_CB", "Phosphates_CB", "CO3_over_Phosphates_CB",
                            "Phosphates_over_OM_MB", "Phosphates_MB", "CO3_over_Phosphates_MB",
                            "PC1", "PC2", "PC3")

combined_tga_ct_fits_main_phenotypes <- filter(combined_tga_ct_fits,
                                               variable %in% main_tga_ct_phenotypes)

subplot_order <- match(main_tga_ct_phenotypes, unique(combined_tga_ct_fits_main_phenotypes$variable))

combined_tga_ct_fits_main_phenotypes$pretty_name <- 
    gsub(combined_tga_ct_fits_main_phenotypes$pretty_name,
         pattern = "QCT PC1 'high density, thickness, content'",
         replacement = "QCT PC1 'high density,\nthickness, content'")

combined_tga_ct_fits_main_phenotypes$pretty_name <- 
    factor(combined_tga_ct_fits_main_phenotypes$pretty_name,
           levels = unique(combined_tga_ct_fits_main_phenotypes$pretty_name)[subplot_order])

plot_tga_estimates_main_phenotypes <- ggplot() +
    geom_pointrange(aes(colour = breed,
                        y = fit,
                        ymin = lwr,
                        ymax = upr,
                        x = cage.pen),
                    size = 0.7,
                    position = position_dodge(0.5),
                    data = combined_tga_ct_fits_main_phenotypes) +
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




combined_tga_ct_contr <- rbind(contr_ct,
                               tga_contr_pretty)

combined_tga_ct_contr$pretty_name <- factor(combined_tga_ct_contr$pretty_name,
                                            levels = c(pretty_tga_trait_names$pretty_name[c(8, 1, 9, 2, 10, 3, 11,
                                                                                            4, 12, 5, 13, 6, 14, 7)],
                                                       unique(contr_ct$pretty_name)))


plot_tga_contr <- qplot(x = contrast,
                        y = estimate,
                        ymin = conf.low,
                        ymax = conf.high,
                        geom = "pointrange",
                        size = I(0.7),
                        data = combined_tga_ct_contr) +
    facet_wrap(~ pretty_name, scale = "free", ncol = 2) +
    geom_hline(yintercept = 0,
               colour = "red",
               linetype = 2) +
    coord_flip() +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.background = element_blank()) +
    xlab("") +
    ylab("") +
    ggtitle("Difference between housing systems (cage minus pen)")


pdf("figures/plot_tga_ct_coef.pdf", width = 8, height = 10)
print(plot_tga_estimates)
dev.off()

pdf("figures/plot_tga_ct_coef_main_phenotypes.pdf", width = 8, height = 10)
print(plot_tga_estimates_main_phenotypes)
dev.off()

pdf("figures/plot_tga_ct_contr.pdf", width = 8, height = 10)
print(plot_tga_contr)
dev.off()



## Plot of main TGA, CT traits with medullary and cortical on the same panel

subplot_trait <- function(data) {
    ggplot() +
        geom_pointrange(aes(colour = breed,
                            y = fit,
                            ymin = lwr,
                            ymax = upr,
                            x = cage.pen),
                        position = position_dodge(0.5),
                        data = data) +
        facet_wrap(~ pretty_name) +
        scale_colour_manual(values = c("blue", "red"),
                            name = "") +
        theme_bw() +
        theme(panel.grid = element_blank(),
              strip.background = element_blank()) +
        xlab("") +
        ylab("")
}

subplot_trait(filter(combined_tga_ct_fits_main_phenotypes,
                     variable %in% c("CO3_over_Phosphates_CB", "CO3_over_Phosphates_MB")))



## Combined correlation heatmap

pheno_split <- split(pheno[,c(2, 7, ct_ix, tga_ix)], pheno$cage.pen)
cor_split <- lapply(pheno_split, cor, use = "p")

cor_split_df <- lapply(cor_split, function(cors) {
    cors <- data.frame(cors)
    cors$variable1 <- rownames(cors)
    cors
})
cor_split_long <- map_dfr(cor_split_df,
                          pivot_longer,
                          cols = -variable1,
                          .id = "cage.pen",
                          names_to = "variable2")


pretty_trait_names <- rbind(pretty_ct_trait_names,
                            pretty_tga_trait_names,
                            data.frame(name = c("weight", "load_N"),
                                       pretty_name = c("body weight", "bone breaking strength"),
                                       stringsAsFactors = FALSE))

cor_split_long <- inner_join(cor_split_long,
                             pretty_trait_names,
                             by = c("variable1" = "name"))

cor_split_long <- inner_join(cor_split_long,
                             pretty_trait_names,
                             by = c("variable2" = "name"))

colnames(cor_split_long)[5:6] <- c("variable1_pretty", "variable2_pretty")

cor_split_long$variable1_pretty <- factor(cor_split_long$variable1_pretty,
                                         levels = pretty_trait_names$pretty_name)

cor_split_long$variable2_pretty <- factor(cor_split_long$variable2_pretty,
                                          levels = pretty_trait_names$pretty_name)



plot_pqct_tga_split_heatmap <- ggplot() +
    geom_tile(aes(x = variable1_pretty,
                  y = variable2_pretty,
                  fill = value),
              data = cor_split_long) +
    facet_wrap(~ cage.pen, ncol = 1) +
    scale_fill_gradient2(limits = c(-1, 1),
                         name = "Correlation") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          panel.grid = element_blank(),
          strip.background = element_blank()) +
    ggtitle("Correlations between bone phenotypes") +
    xlab("") +
    ylab("")


pdf("figures/plot_bone_heatmap.pdf", width = 8, height = 10)
print(plot_pqct_tga_split_heatmap)
dev.off()



## Random forest classification

library(randomForest)

set.seed(20200925)

rf_pheno <- na.exclude(pheno[, c(3, ct_ix, tga_ix)])

rf_pheno_split <- split(rf_pheno, rf_pheno$cage.pen)

n_cage <- nrow(rf_pheno_split[[1]])
n_pen <- nrow(rf_pheno_split[[2]])

test_cage_ix <- sample(1:n_cage, 50)
test_pen_ix <- sample(1:n_pen, 50)

training_cage <- rf_pheno_split[[1]][-test_cage_ix,]
training_pen <- rf_pheno_split[[2]][-test_pen_ix,]

training <- rbind(training_cage,
                  training_pen)

testing_cage <- rf_pheno_split[[1]][test_cage_ix,]
testing_pen <- rf_pheno_split[[2]][test_pen_ix,]

testing <- rbind(testing_cage,
                 testing_pen)


rfs <- lapply(c(100, 500, 1000),
              function(ntree) 
                  randomForest(cage.pen ~ .,
                               data = training,
                               ntree = ntree,
                               importance = TRUE))


rf_pheno_ct <- na.exclude(pheno[, c(3, ct_ix)])

rf_pheno_ct_split <- split(rf_pheno_ct, rf_pheno_ct$cage.pen)

n_cage_ct <- nrow(rf_pheno_ct_split[[1]])
n_pen_ct <- nrow(rf_pheno_ct_split[[2]])

test_cage_ix_ct <- sample(1:n_cage_ct, 50)
test_pen_ix_ct <- sample(1:n_pen_ct, 50)

training_cage_ct <- rf_pheno_ct_split[[1]][-test_cage_ix_ct,]
training_pen_ct <- rf_pheno_ct_split[[2]][-test_pen_ix_ct,]

training_ct <- rbind(training_cage_ct,
                     training_pen_ct)

testing_cage_ct <- rf_pheno_ct_split[[1]][test_cage_ix_ct,]
testing_pen_ct <- rf_pheno_ct_split[[2]][test_pen_ix_ct,]

testing_ct <- rbind(testing_cage_ct,
                    testing_pen_ct)


rfs_ct <- lapply(c(100, 500, 1000),
                 function(ntree) 
                     randomForest(cage.pen ~ .,
                                  data = training_ct,
                                  ntree = ntree,
                                  importance = TRUE))

get_oob_error <- function(rf) rf$err.rate[nrow(rf$err.rate), 1]

unlist(lapply(rfs, get_oob_error))
unlist(lapply(rfs_ct, get_oob_error))

randomForest::varImpPlot(rfs[[2]])
randomForest::varImpPlot(rfs_ct[[2]])


rf_pred <- predict(rfs[[2]], testing)
rf_pred_table <- table(rf_pred, testing$cage.pen)


rf_pred_ct <- predict(rfs_ct[[2]], testing_ct)
rf_pred_ct_table <- table(rf_pred_ct, testing_ct$cage.pen)
