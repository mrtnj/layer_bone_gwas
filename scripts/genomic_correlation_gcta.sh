#!/bin/bash

PLINK_PATH=~/tools/plink/
GCTA_PATH=~/tools/gcta/gcta_1.93.0beta_mac/bin/

set -eu

## Prepare files for GCTA 

$PLINK_PATH/plink --allow-extra-chr \
		  --chr-set 40 \
		  --file gwas/all \
		  --make-bed \
		  --out genomic_correlation/all \
		  --maf 0.0001 \
		  --chr 1-28

$GCTA_PATH/gcta64 --bfile genomic_correlation/all \
		  --make-grm-bin \
		  --out genomic_correlation/gcta \
		  --autosome-num 28 --autosome


$GCTA_PATH/gcta64 --reml-bivar \
		  --grm-bin genomic_correlation/gcta \
		  --pheno genomic_correlation/pheno.txt \
		  --qcovar genomic_correlation/qc.txt \
		  --out genomic_correlation/gcta_bivar
