#!/bin/bash

## Run GWAS

GEMMA_PATH=~/tools/gemma/

set -eu


source scripts/run_gwas.sh



run_gwas all_residual_weight \
         all \
	     all_residual_weight.fam \
	     NULL