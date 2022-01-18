#!/bin/bash

PLINK_PATH=~/tools/plink/

set -eu


## Make 012 coded genotype text files for GWAS with hglm

$PLINK_PATH/plink --allow-extra-chr \
		  --chr-set 40 \
		  --file gwas/all \
		  --out gwas/all \
		  --maf 0.0001 \
		  --recodeA

$PLINK_PATH/plink --allow-extra-chr \
		  --chr-set 40 \
		  --file gwas/pen \
		  --out gwas/pen \
		  --maf 0.0001 \
		  --recodeA

$PLINK_PATH/plink --allow-extra-chr \
		  --chr-set 40 \
		  --file gwas/cage \
		  --out gwas/cage \
		  --maf 0.0001 \
		  --recodeA
		  
		  
## Crosses separately
		  
$PLINK_PATH/plink --allow-extra-chr \
		  --chr-set 40 \
		  --file gwas/all_bovans \
		  --out gwas/all_bovans \
		  --maf 0.0001 \
		  --recodeA

$PLINK_PATH/plink --allow-extra-chr \
		  --chr-set 40 \
		  --file gwas/cage_bovans \
		  --out gwas/cage_bovans \
		  --maf 0.0001 \
		  --recodeA
		  
$PLINK_PATH/plink --allow-extra-chr \
		  --chr-set 40 \
		  --file gwas/pen_bovans \
		  --out gwas/pen_bovans \
		  --maf 0.0001 \
		  --recodeA


$PLINK_PATH/plink --allow-extra-chr \
		  --chr-set 40 \
		  --file gwas/all_lsl \
		  --out gwas/all_lsl \
		  --maf 0.0001 \
		  --recodeA

$PLINK_PATH/plink --allow-extra-chr \
		  --chr-set 40 \
		  --file gwas/cage_lsl \
		  --out gwas/cage_lsl \
		  --maf 0.0001 \
		  --recodeA
		  
$PLINK_PATH/plink --allow-extra-chr \
		  --chr-set 40 \
		  --file gwas/pen_lsl \
		  --out gwas/pen_lsl \
		  --maf 0.0001 \
		  --recodeA
