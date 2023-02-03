
## Make a table of sample sizes

library(assertthat)
library(dplyr)
library(purrr)
library(readr)
library(tibble)

source("R/hglm_helper_functions.R")

source("R/hglm_gwas_prepare_data_crosses_separate.R")


get_n <- function(group_trait_list) {
    
    map_dfr(group_trait_list,
            function(group_trait) {
                tibble(n = sum(!group_trait$missing))
            },
            .id = "trait")
    
}

scans <- list(pen_bovans, cage_bovans, pen_lsl, cage_lsl)
names(scans) <- c("Bovans, PEN", "Bovans, CAGE", "LSL, PEN", "LSL, CAGE")

n <- map_dfr(scans,
             get_n,
             .id = "scan")

write.csv(n,
          file = "tables/supplementary_table_sample_size.csv",
          quote = FALSE,
          row.names = FALSE)
