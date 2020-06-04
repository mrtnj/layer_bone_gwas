#!/bin/bash

PLINK_PATH=~/tools/plink/
GEMMA_PATH=~/tools/gemma/

set -eu


cd simulations/shared


for REP in {1..50}; do

    ## Convert ped to bed with plink

    $PLINK_PATH/plink \
		  --file shared${REP}_joint \
		  --make-bed \
		  --out shared${REP}_joint \
		  --maf 0.0001
		  
    $PLINK_PATH/plink \
		  --file shared${REP}_e1 \
		  --make-bed \
		  --out shared${REP}_e1 \
		  --maf 0.0001
		  
    $PLINK_PATH/plink \
		  --file shared${REP}_e2 \
		  --make-bed \
		  --out shared${REP}_e2 \
		  --maf 0.0001


    ## Make GRM

    $GEMMA_PATH/gemma -bfile shared${REP}_joint \
		  -gk 2 \
		  -o shared${REP}_joint
		  
    $GEMMA_PATH/gemma -bfile shared${REP}_e1 \
    	  -gk 2 \
		  -o shared${REP}_e1
		  
    $GEMMA_PATH/gemma -bfile shared${REP}_e2 \
		  -gk 2 \
		  -o shared${REP}_e2


    ## GWAS

    $GEMMA_PATH/gemma -bfile shared${REP}_joint \
			  -k output/shared${REP}_joint.sXX.txt \
			  -c shared${REP}_covar_joint.txt \
			  -lmm 4 \
			  -o shared${REP}_joint
			  
    $GEMMA_PATH/gemma -bfile shared${REP}_e1 \
			  -k output/shared${REP}_e1.sXX.txt \
			  -c shared${REP}_covar_e1.txt \
			  -lmm 4 \
			  -o shared${REP}_e1
			  
    $GEMMA_PATH/gemma -bfile shared${REP}_e2 \
			  -k output/shared${REP}_e2.sXX.txt \
			  -c shared${REP}_covar_e2.txt \
			  -lmm 4 \
			  -o shared${REP}_e2
			  
done