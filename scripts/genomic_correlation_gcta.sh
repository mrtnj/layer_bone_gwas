#!/bin/bash

PLINK_PATH=~/tools/plink/
GCTA_PATH=~/tools/gcta/gcta_1.93.0beta_mac/bin/

set -eu

## Load

$PLINK_PATH/plink --allow-extra-chr \
		  --chr-set 40 \
		  --file gwas/all \
		  --make-bed \
		  --out genomic_correlation/load_N/all \
		  --maf 0.0001 \
		  --chr 1-28

$GCTA_PATH/gcta64 --bfile genomic_correlation/all \
		  --make-grm-bin \
		  --out genomic_correlation/load_N/gcta \
		  --autosome-num 28 --autosome

$GCTA_PATH/gcta64 --reml \
		  --grm-bin genomic_correlation/load_N/gcta \
		  --pheno genomic_correlation/load_N/pheno.txt \
		  --qcovar genomic_correlation/load_N/qc.txt \
		  --covar genomic_correlation/load_N/cc.txt \
		  --out genomic_correlation/load_N/gcta_cage

$GCTA_PATH/gcta64 --reml \
		  --grm-bin genomic_correlation/load_N/gcta \
		  --pheno genomic_correlation/load_N/pheno.txt \
		  --mpheno 2 \
		  --qcovar genomic_correlation/load_N/qc.txt \
		  --covar genomic_correlation/load_N/cc.txt \
		  --out genomic_correlation/load_N/gcta_pen

$GCTA_PATH/gcta64 --reml-bivar \
		  --grm-bin genomic_correlation/load_N/gcta \
		  --pheno genomic_correlation/load_N/pheno.txt \
		  --qcovar genomic_correlation/load_N/qc.txt \
		  --covar genomic_correlation/load_N/cc.txt \
		  --out genomic_correlation/load_N/gcta_bivar


## Weight

$PLINK_PATH/plink --allow-extra-chr \
		  --chr-set 40 \
		  --file gwas/all \
		  --make-bed \
		  --out genomic_correlation/weight/all \
		  --maf 0.0001 \
		  --chr 1-28

$GCTA_PATH/gcta64 --bfile genomic_correlation/all \
		  --make-grm-bin \
		  --out genomic_correlation/weight/gcta \
		  --autosome-num 28 --autosome

$GCTA_PATH/gcta64 --reml \
		  --grm-bin genomic_correlation/weight/gcta \
		  --pheno genomic_correlation/weight/pheno.txt \
		  --covar genomic_correlation/weight/cc.txt \
		  --out genomic_correlation/weight/gcta_cage

$GCTA_PATH/gcta64 --reml \
		  --grm-bin genomic_correlation/weight/gcta \
		  --pheno genomic_correlation/weight/pheno.txt \
		  --mpheno 2 \
		  --covar genomic_correlation/weight/cc.txt \
		  --out genomic_correlation/weight/gcta_pen

$GCTA_PATH/gcta64 --reml-bivar \
		  --grm-bin genomic_correlation/weight/gcta \
		  --pheno genomic_correlation/weight/pheno.txt \
		  --covar genomic_correlation/weight/cc.txt \
		  --out genomic_correlation/weight/gcta_bivar
