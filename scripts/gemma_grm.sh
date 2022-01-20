#!/bin/bash

## Make GRM for GEMMA

set -eu


cd gwas


gemma -bfile all \
		  -gk 2 \
		  -o all_grm

## System-separated

gemma -bfile pen \
		  -gk 2 \
		  -o pen_grm

gemma -bfile cage \
		  -gk 2 \
		  -o cage_grm

## Breed-separated

gemma -bfile all_bovans \
		  -gk 2 \
		  -o all_bovans_grm
		  
gemma -bfile all_lsl \
		  -gk 2 \
		  -o all_lsl_grm


gemma -bfile bovans_cage \
		  -gk 2 \
		  -o bovans_cage_grm
		  
gemma -bfile lsl_cage \
		  -gk 2 \
		  -o lsl_cage_grm
		  
		  
gemma -bfile bovans_pen \
		  -gk 2 \
		  -o bovans_pen_grm
		  
gemma -bfile lsl_pen \
		  -gk 2 \
		  -o lsl_pen_grm