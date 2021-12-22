#!/bin/bash

PLINK_PATH=~/tools/plink/

set -eu

## Convert ped to bed with plink

cd gwas

ln -s map.map all.map
ln -s map.map cage.map
ln -s map.map pen.map

ln -s map.map all_bovans.map
ln -s map.map all_lsl.map

ln -s map.map pen_bovans.map
ln -s map.map pen_lsl.map

ln -s map.map cage_bovans.map
ln -s map.map cage_lsl.map

cd ..


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

## Cross separated

$PLINK_PATH/plink --allow-extra-chr \
		  --chr-set 40 \
		  --file gwas/all_bovans \
		  --make-bed \
		  --out gwas/all_bovans \
		  --maf 0.0001	
		  
$PLINK_PATH/plink --allow-extra-chr \
		  --chr-set 40 \
		  --file gwas/all_lsl \
		  --make-bed \
		  --out gwas/all_lsl \
		  --maf 0.0001	

$PLINK_PATH/plink --allow-extra-chr \
		  --chr-set 40 \
		  --file gwas/pen_bovans \
		  --make-bed \
		  --out gwas/pen_bovans \
		  --maf 0.0001		  

$PLINK_PATH/plink --allow-extra-chr \
		  --chr-set 40 \
		  --file gwas/cage_bovans \
		  --make-bed \
		  --out gwas/cage_bovans \
		  --maf 0.0001		  

$PLINK_PATH/plink --allow-extra-chr \
		  --chr-set 40 \
		  --file gwas/pen_LSL \
		  --make-bed \
		  --out gwas/pen_LSL \
		  --maf 0.0001		  

$PLINK_PATH/plink --allow-extra-chr \
		  --chr-set 40 \
		  --file gwas/cage_LSL \
		  --make-bed \
		  --out gwas/cage_LSL \
		  --maf 0.0001		  

