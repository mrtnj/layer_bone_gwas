
## Simulate a smilar GWAS population

library(AlphaSimR)


## Create a Macs command, based on GENERIC but with more splits

nInd <- 80

split <- 100

genLen <- 1

Ne <- 100

speciesParams <- "1E8 -t 1E-5 -r 4E-6"

speciesHist <- "-eN 0.25 5.0 -eN 2.50 15.0 -eN 25.00 60.0 -eN 250.00 120.0 -eN 2500.00 1000.0"

popSize <- 2 * nInd

splitI <- paste(" -I 4", popSize/4, popSize/4, popSize/4, popSize/4)
splitJ1 <- paste(" -ej", split/(4 * Ne) + 1e-06, "2 1")
splitJ2 <- paste(" -ej", split/(4 * Ne) + 2e-06, "4 3")
splitJ3 <- paste(" -ej", split/(4 * Ne) + 3e-06, "3 1")

macs_command <- paste0(speciesParams, splitI, 
                       " ", speciesHist, splitJ1, splitJ2, splitJ3)

founderpop <- runMacs(nInd = 80,
                      nChr = 10,
                      segSites = 1100,
                      manualCommand = macs_command,
                      manualGenLen = genLen)



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


line1 <- founders[1:20]
line2 <- founders[21:40]
line3 <- founders[41:60]
line4 <- founders[61:80]


## Make crossbreds

crossbred1 <- randCross2(line1,
                         line2,
                         nCrosses = 40,
                         nProgeny = 10,
                         simParam = simparam)

crossbred2 <- randCross2(line3,
                         line4,
                         nCrosses = 40,
                         nProgeny = 10,
                         simParam = simparam)


crossbreds <- c(crossbred1,
                crossbred2)


pca_crossbreds <- prcomp(pullSegSiteGeno(crossbreds, simParam = simparam))
plot(pca_crossbreds$x[,1], pca_crossbreds$x[,2])


## Export genotypes and phenotypes

geno <- pullSnpGeno(crossbreds,
                    simParam = simparam)


writePlink(crossbreds,
           "simulations/simulation1",
           simParam = simparam)


covar <- data.frame(intercept = 1,
                    crossbred = c(rep(0, 400), rep(1, 400)))

write.table(covar,
            file = "simulations/simulation1_covar.txt",
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
            file = "simulations/simulation1_qtl.txt",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)


