#!/bin/bash

PLINK_PATH=~/tools/plink/

set -eu

## Convert ped to bed with plink

ln -s map.map gwas/all.map
ln -s map.map gwas/cage.map
ln -s map.map gwas/pen.map
ln -s map.map gwas/bovans_pen.map
ln -s map.map gwas/bovans_cage.map
ln -s map.map gwas/lsl_pen.map
ln -s map.map gwas/lsl_cage.map


$PLINK_PATH/plink --allow-extra-chr \
		  --chr-set 40 \
		  --file gwas/all \
		  --make-bed \
		  --out gwas/all \
		  --maf 0.0001

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

$PLINK_PATH/plink --allow-extra-chr \
		  --chr-set 40 \
		  --file gwas/bovans_cage \
		  --make-bed \
		  --out gwas/bovans_cage \
		  --maf 0.0001

$PLINK_PATH/plink --allow-extra-chr \
		  --chr-set 40 \
		  --file gwas/bovans_pen \
		  --make-bed \
		  --out gwas/bovans_pen \
		  --maf 0.0001

$PLINK_PATH/plink --allow-extra-chr \
		  --chr-set 40 \
		  --file gwas/lsl_cage \
		  --make-bed \
		  --out gwas/lsl_cage \
		  --maf 0.0001

$PLINK_PATH/plink --allow-extra-chr \
		  --chr-set 40 \
		  --file gwas/lsl_pen \
		  --make-bed \
		  --out gwas/lsl_pen \
		  --maf 0.0001
