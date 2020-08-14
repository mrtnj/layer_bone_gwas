#!/bin/bash

PLINK_PATH=~/tools/plink/

set -eu


## Make 012 coded genotype text files

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
