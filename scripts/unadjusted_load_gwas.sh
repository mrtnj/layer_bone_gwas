#!/bin/bash

## Run GWAS

GEMMA_PATH=~/tools/gemma/

set -eu


source scripts/run_gwas.sh


for GROUP in all pen cage; do

    run_gwas ${GROUP}_unadjusted_load \
             ${GROUP} \
	         fam_${GROUP}_load_N.fam \
	         covar_${GROUP}_breed.txt
	        
done
