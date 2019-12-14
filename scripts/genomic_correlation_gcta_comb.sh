#!/bin/bash

GCTA_PATH=~/tools/gcta/gcta_1.93.0beta_mac/bin/

set -eu



$GCTA_PATH/gcta64 --reml \
		  --grm-bin genomic_correlation/gcta \
		  --pheno genomic_correlation/pheno_comb.txt \
		  --qcovar genomic_correlation/qc.txt \
		  --covar genomic_correlation/cc.txt \
		  --out genomic_correlation/gcta_comb_cage

$GCTA_PATH/gcta64 --reml \
		  --grm-bin genomic_correlation/gcta \
		  --pheno genomic_correlation/pheno_comb.txt \
		  --mpheno 2 \
		  --qcovar genomic_correlation/qc.txt \
		  --covar genomic_correlation/cc.txt \
		  --out genomic_correlation/gcta_comb_pen

$GCTA_PATH/gcta64 --reml-bivar \
		  --grm-bin genomic_correlation/gcta \
		  --pheno genomic_correlation/pheno_comb.txt \
		  --qcovar genomic_correlation/qc.txt \
		  --covar genomic_correlation/cc.txt \
		  --out genomic_correlation/gcta_comb_bivar
