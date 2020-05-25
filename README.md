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

* R/plot_phenotypes.R -- Plots and linear models of phentotypes


## Format plink files for GEMMA

* R/make_pheno_table.R -- Collect all phenotypes in one file

* R/prepare_plink_files.R -- Reads text files of genotypes and phenotypes
to create ped and phenotype files for GEMMA.

* scripsts/convert_bed.sh -- Use plink to convert binary files for GEMMA.


## Run GEMMA for GWAS

* scripts/gemma_grm.sh -- Create genomic relationship matrix

* scripts/gemma_gwas.sh -- Run GWAS

* R/gwas_summary.R -- Summarise and visualise results


## Run GCTA for genomic correlations

* scripts/genomic_correlation_gcta.sh