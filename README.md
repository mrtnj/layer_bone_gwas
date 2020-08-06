# Layer bone GWAS

GWAS of the Lövstad 2015 layer experiment.


## Update SNP chip markers to Galgal6

* R/prepare_marker_fasta.R -- Read SNP chip design file and output source
sequences for alignment.

* scripts/markers_run_blat.sh -- Use Blat to align sequences to Galgal6

* R/markers_filter_alignment.R -- Read output from Blat mand make a map file with
new positions.


## Exploratory analysis

* R/plot_genotypes.R -- QC of SNP genotypes

* R/plot_phenotypes.R -- Plots and linear models of bone strenght and body weight

* R/plot_pqct_tga_phenotypes.R -- Plots and linear models of pQCT and TGA phenotypes

* R/trait_modelling_functions.R -- Helper functions


## Format plink files for GEMMA

* R/make_pheno_table.R -- Collect all phenotypes in one file

* R/prepare_plink_files.R -- Reads text files of genotypes and phenotypes
to create ped and phenotype files for GEMMA.

* scripsts/convert_bed.sh -- Use plink to convert binary files for GEMMA.


## Run GEMMA for GWAS

* scripts/gemma_grm.sh -- Create genomic relationship matrix

* scripts/gemma_gwas.sh -- Run GWAS


## GWAS results

* R/gwas_summary.R -- Summarise and visualise results

* R/gwas_overlap.R -- Investigate overlaps between GWAS for different traits and genes


## Conditional GWAS on the body mass loci on chr4

* R/prepare_conditional_gwas.R -- Make covariate files for conditional GWAS

* scripts/conditional_gwas.sh -- Run conditional GWAS

* R/conditional_gwas_summary.R -- Summarise results of conditional GWAS


## Soundness check of body weight covariate for GWAS

* script/unadjusted_load_gwas.sh -- Run the breaking strength GWAS without body weight covarite, to check if this does indeed detect he body weight loci

* R/unadjusted_load_gwas_summary.R -- Summarise results of soundness check


## Run GCTA for genomic correlations

* scripts/genomic_correlation_gcta.sh -- Run bivariate genomic model with GCTA


## Simulations

* R/simulation_shared.R -- Fake data simulation for shared genetic architecture

* R/simulation_gxe.R -- Fake data simulation for GxE

* scripts/simulated_gwas.sh -- Run GWAS on simulated data

* scripts/simulated_gwas_gxe.sh -- Run GWAS on simulated data with GxE

* R/simulation_results.R -- Summarise results
