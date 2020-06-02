
## Simulate a situation with GxE

library(AlphaSimR)
library(dplyr)


source("R/simulation_helper_functions.R")

founderpop <- create_founder_population()



## Set up genetic architecture

simparam <- SimParam$new(founderpop)

simparam$setGender("no")

simparam$addTraitA(nQtlPerChr = 5,
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


pca_founders <- prcomp(pullSegSiteGeno(founders, simParam = simparam))
plot(pca_founders$x[,1], pca_founders$x[,2])


crossbreds <- make_crossbreds(founders)

pca_crossbreds <- prcomp(pullSegSiteGeno(crossbreds, simParam = simparam))
plot(pca_crossbreds$x[,1], pca_crossbreds$x[,2])


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
           "simulations/gxe/gxe1_joint",
           trait = 3,
           simParam = simparam)

writePlink(crossbreds[which(covar$environment == 0)],
           "simulations/gxe/gxe1_e1",
           trait = 3,
           simParam = simparam)

writePlink(crossbreds[which(covar$environment == 1)],
           "simulations/gxe/gxe1_e2",
           trait = 3,
           simParam = simparam)


write.table(covar,
            file = "simulations/gxe/gxe1_covar_joint.txt",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

write.table(filter(covar, environment == 0)[,c("intercept", "crossbred")],
            file = "simulations/gxe/gxe1_covar_e1.txt",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

write.table(filter(covar, environment == 1)[,c("intercept", "crossbred")],
            file = "simulations/gxe/gxe1_covar_e2.txt",
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
            file = "simulations/gxe/gxe1_qtl.txt",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)


