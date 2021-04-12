#!/bin/bash

## Make GRM for GEMMA

GEMMA_PATH=~/tools/gemma/

set -eu


cd gwas


$GEMMA_PATH/gemma -bfile all \
		  -gk 2 \
		  -o all_grm

## System-separated

$GEMMA_PATH/gemma -bfile pen \
		  -gk 2 \
		  -o pen_grm

$GEMMA_PATH/gemma -bfile cage \
		  -gk 2 \
		  -o cage_grm

