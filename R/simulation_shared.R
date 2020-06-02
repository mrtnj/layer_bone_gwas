
## Simulate a shared genetic architecture where environment doesn't matter

library(AlphaSimR)
library(dplyr)

source("R/simulation_helper_functions.R")


for (rep_ix in 1:10) {

    founderpop <- create_founder_population()


    ## Set up genetic architecture
    
    simparam <- SimParam$new(founderpop)
    
    simparam$setGender("no")
    
    simparam$addTraitA(nQtlPerChr = 5)
    simparam$setVarE(h2 = 0.3)
    
    simparam$addSnpChip(nSnpPerChr = 1000)
    
    founders <- newPop(founderpop,
                       simParam = simparam)
    
    
    pca_founders <- prcomp(pullSegSiteGeno(founders, simParam = simparam))
    plot(pca_founders$x[,1], pca_founders$x[,2])
    
    
    crossbreds <- make_crossbreds(founders)
    
    
    pca_crossbreds <- prcomp(pullSegSiteGeno(crossbreds, simParam = simparam))
    plot(pca_crossbreds$x[,1], pca_crossbreds$x[,2])
    
    
    ## Set up covariates -- environment makes no difference to trait
    
    covar <- data.frame(intercept = 1,
                        crossbred = c(rep(0, 400), rep(1, 400)),
                        environment = c(rep(c(0, 1),
                                            times = 2,
                                            each = 200)))
    
    
    
    ## Export genotypes and phenotypes
    
    geno <- pullSnpGeno(crossbreds,
                        simParam = simparam)
    
    
    writePlink(crossbreds,
               paste("simulations/shared/shared", rep_ix, "_joint", sep = ""),
               simParam = simparam)
    
    writePlink(crossbreds[which(covar$environment == 0)],
               paste("simulations/shared/shared", rep_ix, "_e1", sep = ""),
               simParam = simparam)
    
    writePlink(crossbreds[which(covar$environment == 1)],
               paste("simulations/shared/shared", rep_ix, "_e2", sep = ""),
               simParam = simparam)
    
    
    
    write.table(covar,
                file = paste("simulations/shared/shared", rep_ix,
                             "_covar_joint.txt", sep = ""),
                row.names = FALSE,
                col.names = FALSE,
                quote = FALSE)
    
    write.table(filter(covar, environment == 0)[,c("intercept", "crossbred")],
                file = paste("simulations/shared/shared", rep_ix,
                             "_covar_e1.txt", sep = ""),
                row.names = FALSE,
                col.names = FALSE,
                quote = FALSE)
    
    write.table(filter(covar, environment == 1)[,c("intercept", "crossbred")],
                file = paste("simulations/shared/shared", rep_ix,
                             "_covar_e2.txt", sep = ""),
                row.names = FALSE,
                col.names = FALSE,
                quote = FALSE)
    
    
    
    ## Get information about QTL
    
    qtl_location <- getQtlMap(1, simParam = simparam)
    
    qtl_location$pos <- round(qtl_location$pos * 10^8)
    
    qtl_geno <- pullQtlGeno(crossbreds, simParam = simparam)
    qtl_location$frequency <- colSums(qtl_geno)/nrow(qtl_geno)/2
    
    qtl_location$a <- simparam$traits[[1]]@addEff
    
    qtl_location$variance_explained <- 2 * qtl_location$a ^ 2 *
        qtl_location$frequency * (1 - qtl_location$frequency)
    
    qtl_location$percent_variance_explained <- 100 *
        qtl_location$variance_explained / as.vector(varP(crossbreds))
    
    
    write.table(qtl_location,
                file = paste("simulations/shared/shared", rep_ix,
                             "_qtl.txt", sep = ""),
                row.names = FALSE,
                col.names = TRUE,
                quote = FALSE)
    
    
}