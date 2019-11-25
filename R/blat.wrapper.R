## A small wrapper around standalone blat.


psl.header <- c("matches", "misMatches", "repMatches", "nCount",
                "qNumInsert", "qBaseInsert", "tNumInsert",
                "tBaseInsert", "strand", "qName", "qSize", "qStart",
                "qEnd", "tName", "tSize", "tStart", "tEnd",
                "blockCount", "blockSizes", "qStarts", "tStarts")

percent.identity.score <- function(psl) {
  ## This function calculates the "milliBad" measure of the UCSC
  ## genome browser. It is essentially mismatches per million bases,
  ## but with an additional gap penalty.

  ## 100 - milliBad * 0.1 gives the percentage ID score of an alignment.
  milli.bad <- numeric(nrow(psl))
  for (j in 1:nrow(psl)) {
    q.alignment.size <- psl$qEnd[j] - psl$qStart[j]
    t.alignment.size <- psl$tEnd[j] - psl$tStart[j]
    alignment.size <- min(q.alignment.size, t.alignment.size)
    size.difference <- q.alignment.size - t.alignment.size
    if (size.difference < 0) {
      ## if target alignment is longer than query alignment
      size.difference <-  0
    }
    insert.factor <- psl$qNumInsert[j]
    ## total number of bases aligned
    total <- psl$match[j] + psl$repMatch[j] + psl$misMatch[j]
    if (total != 0) {
      milli.bad[j] <- (1000 * (psl$misMatch[j] + insert.factor +
                               round(3*log(1+size.difference)))) / total
    }
  }
  return(100 - milli.bad * 0.1)
}

alignment.score <- function(psl) {
  ## The score is simply the number of matches minus mismatches and
  ## the number of inserts (in both query and target).
  return(psl$match - psl$misMatches - psl$qNumInsert - psl$tNumInsert)
}


read.psl <- function(filename) {
  psl <- read.table(filename, header=F, skip=5, stringsAsFactors=F)
  colnames(psl) <- psl.header
  psl$percentID <- percent.identity.score(psl)
  psl$score <- alignment.score(psl)
  return(psl)
}



#### Functions to find good alignments.

best.alignment <- function(psl, score.cutoff=0, id.cutoff=0) {
  
  ## Finds the best alignment for each sequence out of a set of blat
  ## results. Score and percent ID cutoffs can be applied.
  
  n <- nrow(psl)
  best.psl <- data.frame(matches=integer(n),
                         misMatches=integer(n),
                         repMatches=integer(n),
                         nCount=integer(n),
                         qNumInsert=integer(n),
                         qBaseInsert=integer(n),
                         tNumInsert=integer(n),
                         tBaseInsert=integer(n),
                         strand=character(n),
                         qName=character(n),
                         qSize=integer(n),
                         qStart=integer(n),
                         qEnd=integer(n),
                         tName=character(n),
                         tSize=integer(n),
                         tStart=integer(n),
                         tEnd=integer(n),
                         blockCount=integer(n),
                         blockSizes=character(n),
                         qStarts=character(n),
                         tStarts=character(n),
                         percentID=numeric(n),
                         score=integer(n),
                         stringsAsFactors=F)

  ## remove alignments that don't pass thresholds
  below.score.cutoff <- which(psl$score < score.cutoff)
  below.id.cutoff <- which(psl$percentID < id.cutoff)
  to.remove <- union(below.score.cutoff, below.id.cutoff)
  if (length(to.remove) > 0) {
    psl <- psl[-union(below.score.cutoff, below.id.cutoff),]
  }

  ## find the best for each query name
##  k <- 1
##  for (j in unique(psl$qName)) {
##    hits <- subset(psl, qName==j)
##    best.hits <- hits[which(hits$score == max(hits$score)),]
##    best.psl[k:(k+nrow(best.hits)-1),] <- best.hits 
##    k <- k + nrow(best.hits)
##  }

    ##return(best.psl[1:(k-1),])

    dplyr::do(dplyr::group_by(psl, qName), {
        .[which(.$score == max(.$score)),]
    })
}






##### Functions to handle FASTA files.

extract.from.fasta <- function(infile, outfile, keywords) {
  ## Extracts from the input file all sequences that have fasta
  ## headers containing any of the keywords, and writes them to the
  ## output file.

  ## NB! Both files should fit in memory at the same time, or this will
  ## be very inefficient.

  infile <- scan(file=infile, sep="\n", what=character())
  outlist <- list()
  include.sequence <- FALSE

  for (j in 1:length(infile)) {
    if (substr(infile[j], 1, 1) == ">") {
      ## fasta header!
      keyword.matches <- unlist(lapply(keywords, grepl, x=infile[j]))
      if (sum(keyword.matches) > 0) {
        ## if it matches any keyword, start including lines
        include.sequence <- TRUE
        outlist <- append(outlist, infile[j])
      } else {
        ## if it doesn't, stop including lines
        include.sequence <- FALSE
      }
    } else {
      if (include.sequence) {
        outlist <- append(outlist, infile[j])
      }
    }

  }
  write(unlist(outlist), file=outfile, sep="\n")
}



