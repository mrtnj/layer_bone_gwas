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


for CHR in {1..28} 30 31 32 33; do
     plink --file gwas/all \
           --chr $CHR \
           --recode \
           --allow-extra-chr \
           --chr-set 40 \
           --geno 0.1 \
           --out gwas/phasing/input/chr$CHR ;
done 


cd gwas/phasing


for CHR in {1..28} 30 31 32 33; do

    shapeit \
        --input-ped input/chr$CHR.ped input/chr$CHR.map
        --output-max output/chr$CHR.phased.haps output/chr$CHR.phased.sample \
        --input-map ../../annotation/elferink2010_shapeit_GRCg6a/chr$CHR.txt \
        --effective-size 1000 \
        --thread 8
        
done
    



