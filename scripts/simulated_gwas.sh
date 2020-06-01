#!/bin/bash

PLINK_PATH=~/tools/plink/
GEMMA_PATH=~/tools/gemma/

set -eu

## Convert ped to bed with plink

$PLINK_PATH/plink \
		  --file simulations/simulation1 \
		  --make-bed \
		  --out simulations/simulation1 \
		  --maf 0.0001
		  
cd simulations


$GEMMA_PATH/gemma -bfile simulation1 \
		  -gk 2 \
		  -o simulation1


$GEMMA_PATH/gemma -bfile simulation1 \
			  -k output/simulation1.sXX.txt \
			  -c simulation1_covar.txt \
			  -lmm 4 \
			  -o simulation1
			  
