
find_detected_qtl <- function(gwas,
                              threshold,
                              flank,
                              qtl) {
    
    suggestive <- filter(gwas, p_wald < threshold)
    
    suggestive_ranges <- GRanges(seqnames = suggestive$chr,
                                 ranges = IRanges(suggestive$ps - 0.5e6,
                                                  suggestive$ps + 0.5e6))
    
    suggestive_regions <- reduce(suggestive_ranges)
    
    
    qtl_ranges <- GRanges(seqnames = qtl$chr,
                          ranges = IRanges(qtl$pos, qtl$pos))
    
    detected <- rep(FALSE, nrow(qtl))
    detected[queryHits(findOverlaps(qtl_ranges, suggestive_regions))] <- TRUE
    
    detected
}


## Create a Macs command, based on GENERIC but with more splits

create_founder_population <- function() {
    
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
    
    runMacs(nInd = 80,
            nChr = 10,
            segSites = 1100,
            manualCommand = macs_command,
            manualGenLen = genLen)
}


make_crossbreds <- function(founders) {
    
    line1 <- founders[1:20]
    line2 <- founders[21:40]
    line3 <- founders[41:60]
    line4 <- founders[61:80]
    
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
    
    crossbreds
}