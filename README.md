# Layer bone GWAS

GWAS of the Lövstad 2015 layer experiment.


## Update SNP chip markers to Galgal6

* R/prepare_marker_fasta.R -- Read SNP chip design file and output source
sequences for alignment.

* scripts/markers_run_blat.sh -- Use Blat to align sequences to Galgal6

* R/markers_filter_alignment.R -- Read output from Blat and make a map file with
new positions.


## Format plink files

* R/make_pheno_table.R -- Collect all phenotypes in one file

* R/prepare_plink_files.R -- Reads text files of genotypes and phenotypes
to create ped and phenotype files


## Exploratory analysis

* R/plot_genotypes.R -- QC of SNP genotypes

* R/plot_phenotypes.R -- Plots and linear models of bone strength and body weight

* R/plot_pqct_tga_phenotypes.R -- Plots and linear models of pQCT and TGA phenotypes

* R/trait_modelling_functions.R -- Helper functions


## Run quantitative genetics and GWAS with hglm

* scripts/convert_bed.sh -- Convert binary plink files for GEMMA

* script/gemma_grm.sh -- Estimate GRM with GEMMA

* R/hglm_gwas_prepare_data.R -- Set up data for main bone and body weight GWAS

* R/hglm_gwas_prepare_data_crosses_separate.R -- Separate crossbred analysis

* R/hglm_gwas_summary_crosses_separate.R -- Main analysis of GWAS of bone and
body weight with separate crossbreds

* R/hglm_gwas_summary_pqct_crosses_separate.R -- Main analysis of GWAS of QCT 
phenotypes with separate crossbreds

* R/hglm_gwas.R -- Old joint GWAS analysis of bone and body weight

* R/hglm_gwas_tga_prepare_data.R -- Set up data for QCT and TGA GWAS

* R/hglm_gwas_tga.R -- QCT and TGA GWAS

* R/hglm_models.R -- Quantitative genetics models with hglm

* R/hglm_models_tga.R -- Quantitative genetics models for QCT and TGA


## GWAS results

* R/R/hglm_gwas_overlap_crosses_separate.R -- Summarise main GWAS results

* R/hglm_gwas_summary.R -- Old joint GWAS analysis

* R/hglm_gwas_summary_pqct_tga.R -- Summarise QCT and TGA GWAS

* R/hglm_gwas_overlap.R -- Overlap GWAS results between traits

* R/candidate_loci.R -- Look at candidate loci


## Run GCTA for genomic correlations

* R/genomic_correlation_prepare_files.R -- Prepare data for genomic correlation with GCTA

* scripts/genomic_correlation_gcta.sh -- Run bivariate genomic model with GCTA



