#!/bin/bash
#This script adds the PG and PL fields to STITCH-generated output using AWK.
#Usage: 
#./add_PG_PL.sh <STITCH.vcf.gz>

i=$1
tabix $i;
bcftools view $i -h > ${i/.vcf.gz/.header}; 
awk '/^#/;/reference and alternate alleles/ {print "##INFO=<ID=INDEL,Number=0,Type=Flag,Description=\"\">"};/Dosage/ {print "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">";print "##FORMAT=<ID=PG,Number=1,Type=String,Description=\"Best Guessed Genotype with posterior probability threshold of 0.9\">";}' ${i/.vcf.gz/.header}  >  ${i/.vcf.gz/.header}.1; mv ${i/.vcf.gz/.header}.1 ${i/.vcf.gz/.header}

cp ${i/.vcf.gz/.header} ${i/.vcf.gz/.PL.vcf}; bcftools view $i -H | awk 'BEGIN {OFS="\t"};{$9="GT:PG:GP:DS:PL"; for (i=10; i<=NF; i++) {split($i,field,":");pg=field[1];split(field[2],gp,",");if(gp[1] > gp[2] && gp[1]>gp[3]) {field[1]="0/0"};if(gp[2] > gp[1] && gp[2]>gp[3]){field[1]="0/1"};if(gp[3]>gp[1] && gp[3]>gp[2]){field[1]="1/1"};pl1=-log(gp[1])/log(10)*10;pl2=-log(gp[2])/log(10)*10;pl3=-log(gp[3])/log(10)*10;if(pl1=="inf") {pl1=255}; if(pl2=="inf") {pl2=255};if(pl3=="inf") {pl3=255};$i=field[1]":"pg":"field[2]":"field[3]":"int(pl1+0.5)","int(pl2+0.5)","int(pl3+0.5)};print $0}' >> ${i/.vcf.gz/.PL.vcf}; bgzip ${i/.vcf.gz/.PL.vcf}; tabix ${i/.vcf.gz/.PL.vcf.gz};
