## Parse blat output to place SNPs on Galgal6

library(dplyr)

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


## Remove duplicates and poorly aligned SNPs
duplicate_alignments <- new_snps$ID[which(duplicated(new_snps$ID))]


new_snps_pruned <- filter(new_snps,
                          !ID %in% duplicate_alignments &
                          blockCount < 3)



new_snps_pruned$new.pos <- unlist(plyr::alply(new_snps_pruned, 1, find.new.position))
new_snps_pruned$new.chr <- new_snps_pruned$tName

write.csv(new_snps_pruned, file = "outputs/map_Galgal6.csv")


map <- data.frame(chr = new_snps_pruned$new.chr,
                  identifier = new_snps_pruned$ID,
                  genetic_position_dummy = 0,
                  position_bp = new_snps_pruned$new.pos,
                  stringsAsFactors = FALSE)



map_numeric_chr <- filter(map, chr %in% 1:33)
map_numeric_chr <- map_numeric_chr[order(map_numeric_chr$chr,
                                         map_numeric_chr$position_bp),]

map_zwm <- filter(map, chr %in% c("Z", "W", "M"))
map_zwm <- map_zwm[order(map_zwm$chr,
                         map_zwm$position_bp),]

map_un <- filter(map, grepl("Un_", chr))
map_un <- map_un[order(map_un$chr,
                       map_un$position_bp),]

plink_map <- rbind(map_numeric_chr,
                   map_zwm,
                   map_un)

write.table(plink_map,
            file = "gwas/cage.map",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

write.table(plink_map,
            file = "gwas/pen.map",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)
