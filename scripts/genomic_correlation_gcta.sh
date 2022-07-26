#!/bin/bash

PLINK_PATH=~/tools/plink/
GCTA_PATH=~/tools/gcta/

set -eu

## Load, both crosses combined

$PLINK_PATH/plink --allow-extra-chr \
		  --chr-set 40 \
		  --file gwas/all \
		  --make-bed \
		  --out genomic_correlation/load_N/all \
		  --maf 0.0001 \
		  --chr 1-28

$GCTA_PATH/gcta64 --bfile genomic_correlation/load_N/all \
		  --make-grm-bin \
		  --make-grm-alg 1 \
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
		  
		  
## Bovans load

$PLINK_PATH/plink --allow-extra-chr \
		  --chr-set 40 \
		  --file gwas/all_bovans \
		  --make-bed \
		  --out genomic_correlation/load_N_bovans/all_bovans \
		  --maf 0.0001 \
		  --chr 1-28
		  
$GCTA_PATH/gcta64 --bfile genomic_correlation/load_N_bovans/all_bovans \
		  --make-grm-bin \
		  --make-grm-alg 1 \
		  --out genomic_correlation/load_N_bovans/gcta \
		  --autosome-num 28 --autosome
		  
$GCTA_PATH/gcta64 --reml-bivar \
		  --grm-bin genomic_correlation/load_N_bovans/gcta \
		  --pheno genomic_correlation/load_N_bovans/pheno.txt \
		  --qcovar genomic_correlation/load_N_bovans/qc.txt \
		  --out genomic_correlation/load_N_bovans/gcta_bivar
		  
		  
## LSL load

$PLINK_PATH/plink --allow-extra-chr \
		  --chr-set 40 \
		  --file gwas/all_lsl \
		  --make-bed \
		  --out genomic_correlation/load_N_lsl/all_lsl \
		  --maf 0.0001 \
		  --chr 1-28
		  
$GCTA_PATH/gcta64 --bfile genomic_correlation/load_N_lsl/all_lsl \
		  --make-grm-bin \
		  --make-grm-alg 1 \
		  --out genomic_correlation/load_N_lsl/gcta \
		  --autosome-num 28 --autosome
		  
$GCTA_PATH/gcta64 --reml-bivar \
		  --grm-bin genomic_correlation/load_N_lsl/gcta \
		  --pheno genomic_correlation/load_N_lsl/pheno.txt \
		  --qcovar genomic_correlation/load_N_lsl/qc.txt \
		  --out genomic_correlation/load_N_lsl/gcta_bivar
		  


## Weight, both crosses combined

$PLINK_PATH/plink --allow-extra-chr \
		  --chr-set 40 \
		  --file gwas/all \
		  --make-bed \
		  --out genomic_correlation/weight/all \
		  --maf 0.0001 \
		  --chr 1-28

$GCTA_PATH/gcta64 --bfile genomic_correlation/weight/all \
		  --make-grm-bin \
		  --make-grm-alg 1 \
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



## Bovans weight

$PLINK_PATH/plink --allow-extra-chr \
		  --chr-set 40 \
		  --file gwas/all_bovans \
		  --make-bed \
		  --out genomic_correlation/weight_bovans/all_bovans \
		  --maf 0.0001 \
		  --chr 1-28

$GCTA_PATH/gcta64 --bfile genomic_correlation/weight_bovans/all_bovans \
		  --make-grm-bin \
		  --make-grm-alg 1 \
		  --out genomic_correlation/weight_bovans/gcta \
		  --autosome-num 28 --autosome

$GCTA_PATH/gcta64 --reml-bivar \
		  --grm-bin genomic_correlation/weight_bovans/gcta \
		  --pheno genomic_correlation/weight_bovans/pheno.txt \
		  --out genomic_correlation/weight_bovans/gcta_bivar


## LSL weight

$PLINK_PATH/plink --allow-extra-chr \
		  --chr-set 40 \
		  --file gwas/all_lsl \
		  --make-bed \
		  --out genomic_correlation/weight_lsl/all_lsl \
		  --maf 0.0001 \
		  --chr 1-28

$GCTA_PATH/gcta64 --bfile genomic_correlation/weight_bovans/all_lsl \
		  --make-grm-bin \
		  --make-grm-alg 1 \
		  --out genomic_correlation/weight_lsl/gcta \
		  --autosome-num 28 --autosome

$GCTA_PATH/gcta64 --reml-bivar \
		  --grm-bin genomic_correlation/weight_lsl/gcta \
		  --pheno genomic_correlation/weight_lsl/pheno.txt \
		  --out genomic_correlation/weight_lsl/gcta_bivar
