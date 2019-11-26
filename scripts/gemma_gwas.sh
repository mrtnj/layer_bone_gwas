#!/bin/bash

## Run GWAS

GEMMA_PATH=~/tools/gemma/

set -eu



function run_gwas () {
    NAME=$1
    PLINK_BASE=$2
    TRAIT_FILE=$3
    COVAR_FILE=$4

    if [ ! -d gwas/$NAME ]; then
	mkdir gwas/$NAME
    fi

    cd gwas/$NAME

    if [ ! -f $NAME.bed ]; then
	ln -s ../$PLINK_BASE.bim $NAME.bim
	ln -s ../$PLINK_BASE.bed $NAME.bed
	ln -s ../$PLINK_BASE.map $NAME.map
	ln -s ../$TRAIT_FILE $NAME.fam
	if [ COVAR_FILE != "NULL" ]; then
	    ln -s ../$COVAR_FILE covar_$NAME.txt
	fi
    fi

    if [ $COVAR_FILE != "NULL" ]; then
	$GEMMA_PATH/gemma -bfile $NAME \
			  -k ../output/${PLINK_BASE}_grm.sXX.txt \
			  -c covar_$NAME.txt \
			  -lmm 4 \
			  -o $NAME
    fi
    if [ $COVAR_FILE == "NULL" ]; then
	$GEMMA_PATH/gemma -bfile $NAME \
			  -k ../output/${PLINK_BASE}_grm.sXX.txt \
			  -lmm 4 \
			  -o $NAME
    fi
    
    cd ../..
}

run_gwas pen_load_adj \
	 pen \
	 fam_pen_load.fam \
	 covar_pen_breed_weight.txt

run_gwas pen_load \
	 pen \
	 fam_pen_load.fam \
	 covar_pen_breed.txt

run_gwas pen_weight \
	 pen \
	 fam_pen_weight.fam \
	 covar_pen_breed.txt

run_gwas pen_comb_adj \
	 pen \
	 fam_pen_comb.fam \
	 covar_pen_breed_weight.txt


run_gwas cage_load_adj \
	 cage \
	 fam_cage_load.fam \
	 covar_cage_breed_weight.txt

run_gwas cage_load \
	 cage \
	 fam_cage_load.fam \
	 covar_cage_breed.txt

run_gwas cage_weight \
	 cage \
	 fam_cage_weight.fam \
	 covar_cage_breed.txt

run_gwas cage_comb_adj \
	 cage \
	 fam_cage_comb.fam \
	 covar_cage_breed_weight.txt


## Breed-separated

run_gwas bovans_pen_load_adj \
	 bovans_pen \
	 fam_bovans_pen_load.fam \
	 covar_bovans_pen_weight.txt

run_gwas bovans_pen_load \
	 bovans_pen \
	 fam_bovans_pen_load.fam \
	 NULL

run_gwas bovans_cage_load_adj \
	 bovans_cage \
	 fam_bovans_cage_load.fam \
	 covar_bovans_cage_weight.txt

run_gwas bovans_cage_load \
	 bovans_cage \
	 fam_bovans_cage_load.fam \
	 NULL

run_gwas lsl_pen_load_adj \
	 lsl_pen \
	 fam_lsl_pen_load.fam \
	 covar_lsl_pen_weight.txt

run_gwas lsl_pen_load \
	 lsl_pen \
	 fam_lsl_pen_load.fam \
	 NULL

run_gwas lsl_cage_load_adj \
	 lsl_cage \
	 fam_lsl_cage_load.fam \
	 covar_lsl_cage_weight.txt

run_gwas lsl_cage_load \
	 lsl_cage \
	 fam_lsl_cage_load.fam \
	 NULL
