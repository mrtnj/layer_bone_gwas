#!/bin/bash

PLINK_PATH=~/tools/plink/

set -eu

## Convert ped to bed with plink


$PLINK_PATH/plink --allow-extra-chr \
		  --chr-set 40 \
		  --file gwas/pen \
		  --make-bed \
		  --out gwas/pen \
		  --maf 0.0001

$PLINK_PATH/plink --allow-extra-chr \
		  --chr-set 40 \
		  --file gwas/cage \
		  --make-bed \
		  --out gwas/cage \
		  --maf 0.0001
