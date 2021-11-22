
## Plot summaries of haplotypes and allele frequencies

library(assertthat)
library(ggplot2)
library(dplyr)
library(purrr)
library(readr)
library(tibble)



## Read output from Shapeit

chr <- paste("chr", c(1:28, 30:33), sep = "")

haps <- map(paste("gwas/phasing/output/", chr, ".phased.haps", sep = ""),
            read_delim,
            delim = " ",
            col_names = FALSE)


sample_data <- read_delim("gwas/phasing/output/chr1.phased.sample",
                          delim = " ")
sample_data <- sample_data[-1,]


## Add cross label to metadata

pheno <- readRDS("outputs/pheno.Rds")

LSL_ids <- na.exclude(pheno$animal_id[pheno$breed == "LSL"])
bovans_ids <- na.exclude(pheno$animal_id[pheno$breed == "Bovans"])

sample_data$cross <- NA_character_
sample_data$cross[sample_data$ID_1 %in% LSL_ids] <- "LSL"
sample_data$cross[sample_data$ID_1 %in% bovans_ids] <- "Bovans"



## Turn into a data frame of phased observations rather than individuals

observation_data <- sample_data[rep(1:nrow(sample_data), each = 2),]
observation_data$observation_id <- 1:nrow(observation_data)


## Get haplotypes in windows for one phased observation

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



## Get 10-SNP haplotypes 

haps10 <- map(haps,
              make_haplotype_windows_chr,
              window_length = 10)


haps10_LSL <- map(haps10,
                  function(chr) chr[, observation_data$observation_id[observation_data$cross == "LSL"]])

assert_that(all(ncol(haps10_LSL) == sum(observation_data$cross == "LSL")))


haps10_bovans <- map(haps10,
                     function(chr) chr[, observation_data$observation_id[observation_data$cross == "Bovans"]])

assert_that(all(ncol(haps10_bovans) == sum(observation_data$cross == "Bovans")))



## Count the haplotypes within windows on a chromosome

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
        
        counts$freq <- counts$count/sum(counts$count)
        
        window_haplotype_counts[[window_ix]] <- counts
        
        major_haplotype_freq[window_ix] <-
            sum(counts$count[which(counts$count == max(counts$count))]) /
            ncol(hap_strings)
            
    }
    
    list(n_window_haplotypes = n_window_haplotypes,
         window_haplotype_counts = window_haplotype_counts,
         major_haplotype_freq = major_haplotype_freq)
    
}


## Get statistics for the genome

haps10_stats <- map(haps10, 
                    get_window_stats)

haps10_stats_LSL <- map(haps10_LSL,
                        get_window_stats)

haps10_stats_bovans <- map(haps10_bovans,
                           get_window_stats)



## Histogram of the number of haplotypes per window

n_window_haplotypes10 <- unlist(lapply(haps10_stats, "[", "n_window_haplotypes"))
n_window_haplotypes10_LSL <- unlist(lapply(haps10_stats_LSL, "[", "n_window_haplotypes"))
n_window_haplotypes10_bovans <- unlist(lapply(haps10_stats_bovans, "[", "n_window_haplotypes"))

major_haplotype_freq <- unlist(lapply(haps10_stats, "[", "major_haplotype_freq"))
major_haplotype_freq_LSL <- unlist(lapply(haps10_stats_LSL, "[", "major_haplotype_freq"))
major_haplotype_freq_bovans <- unlist(lapply(haps10_stats_bovans, "[", "major_haplotype_freq"))




## Comparison of number of shared haplotypes (Jaccard)

compare_haplotype_sharing <- function(hap_stats_LSL,
                                      hap_stats_bovans) {
    
    n_windows <- length(hap_stats_LSL$window_haplotype_counts)
    
    assert_that(n_windows == length(hap_stats_bovans$window_haplotype_counts))
    
    haplotype_sharing <- numeric(n_windows)
    
    for (window_ix in 1:n_windows) {
        bovans_window <- hap_stats_bovans$window_haplotype_counts[[window_ix]]
        LSL_window <- hap_stats_LSL$window_haplotype_counts[[window_ix]]
        
        n_shared <- length(intersect(bovans_window$haplotype,
                                     LSL_window$haplotype))
        n_total <- length(union(bovans_window$haplotype,
                                LSL_window$haplotype))
        
        haplotype_sharing[window_ix] <- n_shared/n_total
        
    }
    
    haplotype_sharing
}

haplotype_sharing10 <- pmap(list(hap_stats_LSL = haps10_stats_LSL,
                                 hap_stats_bovans = haps10_stats_bovans),
                            compare_haplotype_sharing)



## Comparison of the frequency of each haplotype

compare_haplotype_frequency <- function(hap_stats_LSL,
                                        hap_stats_bovans) {
    
    n_windows <- length(hap_stats_LSL$window_haplotype_counts)
    
    assert_that(n_windows == length(hap_stats_bovans$window_haplotype_counts))
    
    freq_diff <- vector(length = n_windows,
                        mode = "list")
    
    for (window_ix in 1:n_windows) {
        bovans_window <- hap_stats_bovans$window_haplotype_counts[[window_ix]]
        LSL_window <- hap_stats_LSL$window_haplotype_counts[[window_ix]]
        
        counts_combined <- inner_join(bovans_window, ## x
                                      LSL_window, ## y
                                      by = "haplotype")
        colnames(counts_combined)[2:5] <- c("count_bovans", "freq_bovans",
                                            "count_LSL", "freq_LSL")
        
        counts_combined$bovans_minus_LSL <-
            counts_combined$freq_bovans - counts_combined$freq_LSL
        
        freq_diff[[window_ix]] <- counts_combined
    }
    
    freq_diff
}


frequency_difference10 <- pmap(list(hap_stats_LSL = haps10_stats_LSL,
                                    hap_stats_bovans = haps10_stats_bovans),
                               compare_haplotype_frequency)
