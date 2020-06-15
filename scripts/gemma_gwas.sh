#!/bin/bash

## Run GWAS

GEMMA_PATH=~/tools/gemma/

set -eu


source scripts/run_gwas.sh


for TRAIT in load_N comb_g ct_pc1 ct_pc2 ct_pc3 \
             WaterLost_CB OMLost_CB CO2Lost_CB \
             Phosphates_CB Mineral_CB Phosphates_over_OM_CB \
             CO3_over_Phosphates_CB WaterLost_MB OMLost_MB \
             CO2Lost_MB Phosphates_MB Mineral_MB \
             Phosphates_over_OM_MB CO3_over_Phosphates_MB; do
            
             for GROUP in all pen cage; do
            
                 run_gwas ${GROUP}_${TRAIT} \
                          ${GROUP} \
	                      fam_${GROUP}_${TRAIT}.fam \
	                      covar_${GROUP}_breed_weight.txt
	                     
             done
            
done

for GROUP in all pen cage; do

    run_gwas ${GROUP}_weight \
	          ${GROUP} \
	          fam_${GROUP}_weight.fam \
	          covar_${GROUP}_breed.txt
	          
done



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
