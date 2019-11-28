
## Read GWAS results and make graphs

library(dplyr)
library(ggplot2)
library(qqman)
library(readr)
library(stringr)


## Collate results

files <- system("ls gwas/*/output/*assoc.txt", intern = TRUE)

scan_name <- str_match(files, "/([a-z_]+)\\.assoc\\.txt")[,2]

results <- lapply(files,
                  read_tsv,
                  col_types = "ccnnccnnnnnnnnn")

for (file_ix in 1:length(results)) {
    results[[file_ix]]$scan_name <- scan_name[file_ix]
    results[[file_ix]]$numeric_chr <- as.numeric(results[[file_ix]]$chr)
}

## qq(results[[8]]$p_wald)
## manhattan(na.exclude(results[[5]]), chr = "numeric_chr", p = "p_wald", bp = "ps", snp = "rs")



gwas <- Reduce(rbind, results)


saveRDS(gwas,
        file = "outputs/gwas.Rds")



## Plots

chr_lengths <- summarise(group_by(gwas, chr), length = unique(max(ps)))
preferred_order <- c(1:28, 30:31, 33, "Un_NW_020110160v1", "Un_NW_020110165v1", "Z")

chr_lengths <- chr_lengths[match(preferred_order, chr_lengths$chr),]

flatten_coordinates <- function(chr,
                                pos,
                                chr_lengths) {
    pos_flat <- pos
    offset <- 0
 
    for (chr_ix in 1:nrow(chr_lengths)) {
        on_chr <- chr == chr_lengths$chr[chr_ix]
        pos_flat[on_chr] <- pos[on_chr] + offset
        offset <- offset + chr_lengths$length[chr_ix]
    }
 
    pos_flat
}


gwas$global_pos <- flatten_coordinates(gwas$chr, gwas$ps, chr_lengths)

