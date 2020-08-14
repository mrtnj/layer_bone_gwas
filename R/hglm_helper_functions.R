

decompose_grm <- function(grm, missing) {
    svd <- svd(grm[!missing, !missing])
    svd$u %*% diag(sqrt(svd$d))
}


get_covar <- function(pheno,
                      load_data) {
    
    covar <- filter(pheno,
                    animal_id %in% load_data$X1 &
                        !is.na(load_N))[, c("animal_id", "group", "weight", "cage.pen")]
    
    covar <- covar[order(covar$animal_id),]
    
    assert_that(identical(covar$animal_id, load_data$X1))
    
    covar
}


get_h2 <- function(model) {
    model$varRanef[1] / (model$varFix + sum(model$varRanef))
}

get_group_ratio <- function(model) {
    model$varRanef[2] / (model$varFix + sum(model$varRanef))
}