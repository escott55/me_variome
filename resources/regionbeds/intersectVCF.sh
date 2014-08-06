#!/bin/bash

BEDFILE=percentcoverage.bed
VCFFILE=/local/gleeson/resources/onekgenomes/onekgenomes.snps_indels.vcf.gz
#VCFFILE=/home/escott/projects/inbreed/inhouse/allplates/merged_snps.vcf.gz

VCFFILENAME=$(basename $VCFFILE)
EXTENSION=${VCFFILENAME##*.}
VCFFILENAME=${VCFFILENAME%.*}
echo "$VCFFILENAME"

BEDFILENAME=$(basename $BEDFILE)
EXTENSION=${BEDFILENAME##*.}
BEDFILENAME=${BEDFILENAME%.*}
echo "$BEDFILENAME"

echo "$VCFFILENAME.$BEDFILENAME.vcf"
OUTFILE=$VCFFILENAME.$BEDFILENAME.vcf

gunzip -c $VCFFILE | head -300 | awk '$1 ~ /#/' > $OUTFILE
intersectBed -wa -a $VCFFILE -b $BEDFILE >> $OUTFILE

