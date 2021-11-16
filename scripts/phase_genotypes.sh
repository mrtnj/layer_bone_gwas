#!/bin/bash


## Use Shapeit to phase the genotypes


if [ ! -d gwas/phasing ]; then
    mkdir gwas/phasing
fi

if [ ! -d gwas/phasing/input ]; then
    mkdir gwas/phasing/input
fi

if [ ! -d gwas/phasing/output ]; then
    mkdir gwas/phasing/output
fi

## {1..28} 30 31 32 33

for CHR in 28; do
     plink --file gwas/all \
           --chr $CHR \
           --recode \
           -- allow-extra-chr \
           --make-bed \
           --out gwas/phased/input/chr$CHR ;
done 


cd gwas/phasing


for CHR in 28; do

    shapeit \
        --input-bed input/chr$CHR.ped input/chr$CHR.bim input/chr$CHR.fam \
        --output-max output/chr$CHR.phased.haps chr$CHR.phased.sample \
        --input-map ../../annotation/elferink2010_shapeit_GRCg6a/chr$CHR.txt \
        --effective-size 1000 \
        --thread 8
        
done
    



