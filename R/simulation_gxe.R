
## Simulate a situation with GxE

library(AlphaSimR)
library(dplyr)


source("R/simulation_helper_functions.R")


n_rep <- 50

for (rep_ix in 1:n_rep) {

    founderpop <- create_founder_population()


    ## Set up genetic architecture
    
    simparam <- SimParam$new(founderpop)
    
    simparam$setGender("no")
    
    simparam$addTraitA(nQtlPerChr = 1,
                       mean = c(0, 0),
                       var = c(1, 1),
                       corA = matrix(c(1, 0.6,
                                       0.6, 1),
                                     byrow = TRUE,
                                     nrow = 2))
    simparam$setVarE(h2 = c(0.3, 0.3))
    
    simparam$addSnpChip(nSnpPerChr = 1000)
    
    founders <- newPop(founderpop,
                       simParam = simparam)


    crossbreds <- make_crossbreds(founders)
    

    ## Set up covariates
    
    covar <- data.frame(intercept = 1,
                        crossbred = c(rep(0, 400), rep(1, 400)),
                        environment = c(rep(c(0, 1),
                                            times = 2,
                                            each = 200)))
    
    
    ## Fix traits
    
    observed_pheno <- ifelse(covar$environment == 0,
                             crossbreds@pheno[,1],
                             crossbreds@pheno[,2])
    
    
    crossbreds@pheno <- cbind(crossbreds@pheno, observed_pheno)
    
    
    ## Export genotypes and phenotypes
    
    
    writePlink(crossbreds,
               paste("simulations/gxe/gxe", rep_ix, "_joint", sep = ""),
               trait = 3,
               simParam = simparam)
    
    writePlink(crossbreds[which(covar$environment == 0)],
               paste("simulations/gxe/gxe", rep_ix, "_e1", sep = ""),
               trait = 3,
               simParam = simparam)
    
    writePlink(crossbreds[which(covar$environment == 1)],
               paste("simulations/gxe/gxe", rep_ix, "_e2", sep = ""),
               trait = 3,
               simParam = simparam)
    
    
    write.table(covar,
                file = paste("simulations/gxe/gxe", rep_ix, "_covar_joint.txt", sep = ""),
                row.names = FALSE,
                col.names = FALSE,
                quote = FALSE)
    
    write.table(filter(covar, environment == 0)[,c("intercept", "crossbred")],
                file = paste("simulations/gxe/gxe", rep_ix, "_covar_e1.txt", sep = ""),
                row.names = FALSE,
                col.names = FALSE,
                quote = FALSE)
    
    write.table(filter(covar, environment == 1)[,c("intercept", "crossbred")],
                file = paste("simulations/gxe/gxe", rep_ix, "_covar_e2.txt", sep = ""),
                row.names = FALSE,
                col.names = FALSE,
                quote = FALSE)
    
    
    ## Get information about QTL
    
    qtl_location <- getQtlMap(1, simParam = simparam)
    
    qtl_location$pos <- round(qtl_location$pos * 10^8)
    
    qtl_geno <- pullQtlGeno(crossbreds, simParam = simparam)
    qtl_location$frequency <- colSums(qtl_geno)/nrow(qtl_geno)/2
    
    qtl_location$a1 <- simparam$traits[[1]]@addEff
    
    qtl_location$variance_explained1 <- 2 * qtl_location$a1 ^ 2 *
        qtl_location$frequency * (1 - qtl_location$frequency)
    
    qtl_location$percent_variance_explained1 <- 100 *
        qtl_location$variance_explained1 / varP(crossbreds)[1,1]
    
    qtl_location$a2 <- simparam$traits[[2]]@addEff
    
    qtl_location$variance_explained2 <- 2 * qtl_location$a2 ^ 2 *
        qtl_location$frequency * (1 - qtl_location$frequency)
    
    qtl_location$percent_variance_explained2 <- 100 *
        qtl_location$variance_explained2 / varP(crossbreds)[2,2]
    
    
    write.table(qtl_location,
                file = paste("simulations/gxe/gxe", rep_ix, "_qtl.txt", sep = ""),
                row.names = FALSE,
                col.names = TRUE,
                quote = FALSE)
    
    
}