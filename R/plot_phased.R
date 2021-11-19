
## Plot summaries of haplotypes and allele frequencies

library(assertthat)
library(ggplot2)
library(dplyr)
library(purrr)
library(readr)
library(tibble)



chr <- paste("chr", c(1:28, 30:33), sep = "")

haps <- map(paste("gwas/phasing/output/", chr, ".phased.haps", sep = ""),
            read_delim,
            delim = " ",
            col_names = FALSE)





## Get haplotypes in windows for one chromosome

make_haplotype_windows <- function(hap,
                                   window_length) {
    
    n_markers <- length(hap)
    
    window_start <- seq(from = 1,
                        to = n_markers,
                        by = window_length)
    
    n_windows <- length(window_start)
    
    window_end <- window_start + window_length - 1
    window_end[n_windows] <- n_markers
    
    
    hap_strings <- pmap_chr(list(start = window_start,
                                 end = window_end),
                            function(start, end) paste0(hap[start:end],
                                                        collapse = " "))
    
    hap_strings    
}


## Apply to a chromosome

make_haplotype_windows_chr <- function(haps,
                                       window_length) {
    
    haps_no_metadata <- haps[, 6:ncol(haps)]
    
    metadata <- haps[, 1:5]
    
    hap_strings_no_metadata <- map_dfc(haps_no_metadata,
                                       make_haplotype_windows,
                                       window_length = window_length)
    
    
    hap_strings_no_metadata
}



haps10 <- map(haps,
              make_haplotype_windows_chr,
              window_length = 10)


get_window_stats <- function(hap_strings) {
    
    hap_strings <- as.matrix(hap_strings)
    
    n_windows <- nrow(hap_strings)
    
    n_window_haplotypes <- numeric(n_windows)
    major_haplotype_freq <- numeric(n_windows)
    window_haplotype_counts <- vector(length = n_windows,
                                      mode = "list")
    
    for (window_ix in 1:n_windows) {
        
        n_window_haplotypes[window_ix] <- length(unique(hap_strings[window_ix,]))
        
        counts <- as.data.frame(table(hap_strings[window_ix,]))
        colnames(counts) <- c("haplotype", "count")
        
        window_haplotype_counts[[window_ix]] <- counts
        
        major_haplotype_freq[window_ix] <-
            sum(counts$count[which(counts$count == max(counts$count))]) /
            ncol(hap_strings)
            
    }
    
    list(n_window_haplotypes = n_window_haplotypes,
         window_haplotype_counts = window_haplotype_counts,
         major_haplotype_freq = major_haplotype_freq)
    
}


haps10_stats <- map(haps10, 
                    get_window_stats)
