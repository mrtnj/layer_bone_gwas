

## Make a fasta file to move markers to Galgal6

library(Biostrings)
library(readr)



annotation <- read_csv("RA-1698_181010_ResultReport/Chicken_50K_CobbCons_15000986_A.csv",
                       skip = 7)


seq <- annotation$SourceSeq[!is.na(annotation$SourceSeq)]

seq <- sub(seq,
           pattern = "\\[.*\\]",
           replacement = "N")
names(seq) <- annotation$Name[!is.na(annotation$SourceSeq)]

fasta <- DNAStringSet(seq)

writeXStringSet(fasta,
                file = "alignment/source_seq.fasta")
write.table(annotation[, c("Name", "Chr", "MapInfo", "SourceSeq")],
            file = "alignment/source_seq.txt",
            sep = "\t",
            row.names = FALSE)
