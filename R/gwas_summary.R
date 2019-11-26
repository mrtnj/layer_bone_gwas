
## Read GWAS results and make graphs

library(qqman)
library(readr)
library(stringr)


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
