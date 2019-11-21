#!/bin/bash

BLAT_PATH=~/tools/blat/

## Download genome and run blat on F8 SNP sequences

if [ ! -f annotation/galGal6.fa ]; then

    cd annotation
    curl -O https://hgdownload.soe.ucsc.edu/goldenPath/galGal6/bigZips/galGal6.fa.gz
    gunzip galGal6.fa.gz

fi

##$BLAT_PATH/faToTwoBit annotation/galGal6.fa annotation/galGal6.2bit
$BLAT_PATH/blat -t=dna -q=dna annotation/galGal6.2bit alignment/source_seq.fasta alignment/source_seq.psl 
