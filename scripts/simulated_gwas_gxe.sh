#!/bin/bash

PLINK_PATH=~/tools/plink/
GEMMA_PATH=~/tools/gemma/

set -eu

cd simulations/gxe

for REP in {1..10}; do

    ## Convert ped to bed with plink

    $PLINK_PATH/plink \
		  --file gxe1_joint \
		  --make-bed \
		  --out gxe1_joint \
		  --maf 0.0001
		  
    $PLINK_PATH/plink \
		  --file gxe1_e1 \
		  --make-bed \
		  --out gxe1_e1 \
		  --maf 0.0001
		  
    $PLINK_PATH/plink \
		  --file gxe1_e2 \
		  --make-bed \
		  --out gxe1_e2 \
		  --maf 0.0001


    ## Make GRM

    $GEMMA_PATH/gemma -bfile gxe1_joint \
		  -gk 2 \
		  -o gxe1_joint
		  
    $GEMMA_PATH/gemma -bfile gxe1_e1 \
    	  -gk 2 \
		  -o gxe1_e1
		  
    $GEMMA_PATH/gemma -bfile gxe1_e2 \
		  -gk 2 \
		  -o gxe1_e2


    ## GWAS

    $GEMMA_PATH/gemma -bfile gxe1_joint \
			  -k output/gxe1_joint.sXX.txt \
			  -c gxe1_covar_joint.txt \
			  -lmm 4 \
			  -o gxe1_joint
			  
    $GEMMA_PATH/gemma -bfile gxe1_e1 \
			  -k output/gxe1_e1.sXX.txt \
			  -c gxe1_covar_e1.txt \
			  -lmm 4 \
			  -o gxe1_e1
			  
    $GEMMA_PATH/gemma -bfile gxe1_e2 \
			  -k output/gxe1_e2.sXX.txt \
			  -c gxe1_covar_e2.txt \
			  -lmm 4 \
			  -o gxe1_e2
			  
done