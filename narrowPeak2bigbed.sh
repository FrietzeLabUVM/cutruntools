#!/bin/bash
#convert narrowPeak files to bigbed
#assumes hg38
if [ -z $1 ]; then echo need peak file arg1; exit 1; fi
if [ -z $2 ]; then echo need chr sizes file arg2; exit 1; fi
f=$1
f_bed=${f}.2bb.bed
#f_bed=${f_bed/.bed/.2bb.bed}
#f_bed=${f_bed/.narrowPeak/.2bb.bed}
#f_bed=${f_bed/.broadPeak/.2bb.bed}
f_bigbed=${f}.bigbed
#f_bigbed=${f_bigbed/.bed/.bigbed}
#f_bigbed=${f_bigbed/.narrowPeak/.bigbed}
#f_bigbed=${f_bigbed/.broadPeak/.bigbed}
chr_sizes=$2
if [ -f $f_bed ]; then echo bed file $f_bed exists, delete before running.; exit 1; fi
if [ -f $f_bigbed ]; then echo bigbed file $f_bigbed exists, delete before running; exit 1; fi
echo NARROWPEAK is $f
echo BED is $f_bed
echo BIGBED is $f_bigbed
sort -k1,1 -k2,2n $f | awk 'BEGIN {FS="\t"; OFS="\t"} {if ($5 > 1000) $5 = 1000; print $1,$2,$3,$4,$5}' > $f_bed
bedToBigBed $f_bed $chr_sizes $f_bigbed
rm $f_bed
