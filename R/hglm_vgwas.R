library(assertthat)
library(dplyr)
library(hglm)
library(readr)


source("R/hglm_helper_functions.R")


source("R/hglm_gwas_prepare_data.R")

X_disp_null <- model.matrix(~ cage.pen, all_load_covar)

null_model <- hglm(y = all_load$data$X6,
                   X = model.matrix(~ 1 + weight + breed + cage.pen, all_load_covar),
                   Z = cbind(Z_grm_all_load, Z_all_load_group),
                   RandC = c(ncol(Z_grm_all_load), ncol(Z_all_load_group)),
                   X.disp = X_disp_null,
                   calc.like = TRUE)


n_snps <- ncol(snps_pruned_all_load)

marker_id <- character(n_snps)
estimate <- numeric(n_snps)
se <- numeric(n_snps)
LRT <- numeric(n_snps)
p <- numeric(n_snps)

for (snp_ix in 1:n_snps) {
    
    if (snp_ix %% 1000 == 0) {
        cat(snp_ix, "/", n_snps, "\n")
    }

    X_disp <- cbind(model.matrix(~ cage.pen, all_load_covar),
                    snps_pruned_all_load[, snp_ix])

    model <- tryCatch({
        
        hglm(y = all_load$data$X6,
             X = model.matrix(~ 1 + weight + breed + cage.pen, all_load_covar),
             Z = cbind(Z_grm_all_load, Z_all_load_group),
             RandC = c(ncol(Z_grm_all_load), ncol(Z_all_load_group)),
             X.disp = X_disp,
             calc.like = TRUE)

    },
    error = function(condition) {
        message("Error in SNP inex", snp_ix)
        message(condition)
        return(NA)
    })
    
    if (!is.na(model)) {
        marker_id[snp_ix] <- rownames(model$SummVC1)[3]
        estimate[snp_ix] <- model$SummVC1[3, 1]
        se[snp_ix] <- model$SummVC1[3, 2]
        
        test <- lrt(null_model, model)
        LRT[snp_ix] <- test$statistic
        p[snp_ix] <- test$p.value
    } else {
        marker_id[snp_ix] <- colnames(snps_pruned_all_load)[snp_ix]
        estimate[snp_ix] <- NA
        se[snp_ix] <- NA
        LRT[snp_ix] <- NA
        p[snp_ix] <- NA
    }
    
}

vgwas_all_load <- data.frame(marker_id = marker_id,
                             estimate = estimate,
                             se = se,
                             LRT = LRT,
                             p = p,
                             stringsAsFactors = FALSE)


saveRDS(vgwas_all_load,
        file = "gwas/hglm_vgwas_all_load.Rds")