

decompose_grm <- function(grm, missing) {
    svd <- svd(grm[!missing, !missing])
    svd$u %*% diag(sqrt(svd$d))
}


get_covar <- function(pheno,
                      load_data,
                      trait_name) {
    
    covar <- filter(pheno,
                    animal_id %in% load_data$X1 &
                        !is.na( {{ trait_name }} ))[, c("animal_id", "group", "weight", "cage.pen", "breed")]
    
    covar <- covar[match(load_data$X1, covar$animal_id),]
    
    assert_that(identical(covar$animal_id, load_data$X1))
    
    covar
}


get_h2 <- function(model) {
    model$varRanef[1] / (model$varFix + sum(model$varRanef))
}

get_group_ratio <- function(model) {
    model$varRanef[2] / (model$varFix + sum(model$varRanef))
}

get_var_ratio <- function(model, component) {
    model$varRanef[component] / (model$varFix + sum(model$varRanef))
}


## Read plink data

read_data <- function(filename) {
    
    plink_data <- read_delim(filename,
                             col_names = FALSE,
                             delim = " ")
    
    missing <- is.na(plink_data$X6)
    
    list(data = plink_data[!missing,],
         missing = missing)
}


prune_snp_matrix <- function(snp_matrix) {
 
    missing <- lapply(snp_matrix, function(x) sum(is.na(x)))
    
    genotype_values <- lapply(snp_matrix, function(x) length(unique(x)))
    
    snp_matrix[, missing == 0 & genotype_values > 1]
    
}

## Run GWAS

run_gwas <- function(pheno,
                     X,
                     Z,
                     RandC,
                     snp_matrix) {
    
    n_ind <- length(pheno)
    n_snp <- ncol(snp_matrix)
    
    ## Fit baseline model
    
    model_baseline <- hglm(y = pheno,
                           X = X,
                           Z = Z,
                           RandC = RandC)
    
    
    ratio <- model_baseline$varRanef/model_baseline$varFix
    
    V <- RepeatABEL::constructV(Z,
                                RandC,
                                ratio)
    
    eigV <- eigen(V)
    
    transformation_matrix <- diag(1/sqrt(eigV$values)) %*% t(eigV$vectors)
    
    transformed_y <- transformation_matrix %*% pheno
    transformed_X <- transformation_matrix %*% X
    
    
    ## Null model
    
    qr0 <- qr(transformed_X)
    
    est0 <- qr.coef(qr0, transformed_y)
    
    null_residual <- transformed_y - transformed_X %*% est0
    
    RSS_null <- sum(null_residual^2)/n_ind
    
    
    ## SNP models
    
    estimates <- numeric(n_snp)
    LRT <- numeric(n_snp)
    p <- numeric(n_snp)
    
    for (snp_ix in 1:n_snp) {
        
        transformed_snp <- transformation_matrix %*% as.matrix(snp_matrix[, snp_ix])
        
        X1 <- cbind(transformed_snp, transformed_X)
        qr1 <- qr(X1)
        est1 <- qr.coef(qr1, transformed_y)
        residual1 <- transformed_y - X1 %*% est1
        RSS1 <- sum(residual1^2)/n_ind
        
        estimates[snp_ix] <- est1[1]
        LRT[snp_ix] <- -n_ind * (log(RSS1) - log(RSS_null))
        p[snp_ix] <- 1 - pchisq(LRT[snp_ix],
                                df = 1)
    }
    
    data.frame(marker_id = colnames(snp_matrix),
               estimates,
               LRT,
               p)
}