library(broom)
library(multcomp)


model_breed_cagepen <- function(long, bw_covariate = FALSE) {
    
    long_split <- split(long, long$variable)
    
    if (bw_covariate) {
        
        models <- lapply(long_split,
                         function(x) {
                             lm(value ~ breed * cage.pen + weight, data = x)
                         })
        
        levels <- expand.grid(cage.pen = c("CAGE", "PEN"),
                              breed = c("LSL", "Bovans"),
                              weight = mean(long$weight),
                              stringsAsFactors = FALSE)
    } else {
    
        models <- lapply(long_split,
                         function(x) {
                             lm(value ~ breed * cage.pen, data = x)
                     })
        
        levels <- expand.grid(cage.pen = c("CAGE", "PEN"),
                              breed = c("LSL", "Bovans"),
                              stringsAsFactors = FALSE)
    }
    
    
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
    if (bw_covariate) {
    
        cage_minus_pen_within_bovan <- matrix(c(0, 0, -1, mean(long$weight), 0), 1)
        cage_minus_pen_within_lsl <- matrix(c(0, 0, -1, mean(long$weight), -1), 1)
        
    } else {
        
        cage_minus_pen_within_bovan <- matrix(c(0, 0, -1, 0), 1)
        cage_minus_pen_within_lsl <- matrix(c(0, 0, -1, -1), 1)
        
    }
    
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