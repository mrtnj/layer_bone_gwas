## Parse blat output to place SNPs on Galgal6


source("R/blat.wrapper.R")


## Read alignments and select the most high-scoring alignment for each
## sequence.

psl <- read.psl("alignment/source_seq.psl")
nrow(psl)
best_alignments <- best.alignment(psl)
nrow(best_alignments)


## Read the snp list. Find the position of each snp in the given
## sequence.

snps <- read.table("alignment/source_seq.txt", sep = "\t", header = TRUE,
                   stringsAsFactors = FALSE)

find.snp <- function(seq) {
  seq <- strsplit(seq, split="")[[1]]
  which(seq == "[")
}

snps <- na.exclude(snps)

snps$position.in.seq <- unlist(lapply(snps$SourceSeq, find.snp))



## Merge together new and old positions.

new_positions <- best_alignments
new_positions$tName <- gsub(new_positions$tName, pattern = "chr", replacement = "")
new_positions$ID  <- new_positions$qName


new_snps <- merge(new_positions,
                  data.frame(ID = snps$Name,
                             old.chr = snps$Chr, old.pos = snps$MapInfo,
                             position.in.seq = snps$position.in.seq,
                             stringsAsFactors = FALSE),
                  by.x = "ID", by.y = "ID")



## Based on the alignments, find the genomic position of the snp.

## (Assumes each alignment has at most one gap.)

find.new.position <- function(snp) {
  blocks <- as.numeric(strsplit(snp$blockSizes, split=",")[[1]])
  starts <- as.numeric(strsplit(snp$tStarts, split=",")[[1]])
  position <- NA

  if (blocks[1] >= snp$position.in.seq) {
      position <- starts[1] + snp$position.in.seq
  } else if (blocks[1] < snp$position.in.seq) {
      position <- starts[2] + snp$position.in.seq - blocks[1] - 2
  }
  position
}

library(plyr)

new_snps$new.pos <- unlist(alply(new_snps, 1, find.new.position))
new_snps$new.chr <- new_snps$tName

write.csv(new_snps, file = "outputs/map_Galgal6.csv")
