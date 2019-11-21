#!/bin/bash

PLINK_PATH=~/tools/plink/

set -eu

## Convert ped to bed with plink


$PLINK_PATH/plink --file gwas/pen \
		  --make-bed \
		  --out gwas/pen
