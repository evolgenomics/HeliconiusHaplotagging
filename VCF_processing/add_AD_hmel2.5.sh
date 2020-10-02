#!/bin/bash
vcf=$1
chr=`echo $vcf | awk '{split($1,a,".");print substr(a[5],6,6)}'`
#bcftools query -f "%CHROM\t%POS[\t%AD]\n" $vcf -H | bgzip  -c > $vcf.ad.gz
#tabix $vcf.ad.gz -b 2 -e 2 -s 1

echo $vcf"	->	"$chr
#if [[ $vcf  == *1o.mpileup.vcf ]]; then 
bcftools annotate -a stitch.Hmel2.5.Hmel2${chr}.mpileup.vcf.ad.gz  -h ad.hdr -c CHROM,POS,FMT/AD $vcf | bcftools view - -Oz -o  ${vcf/.PL.vcf.gz/.PL.AD.vcf.gz}; #`basename $vcf .PL.vcf.gz`.PL.AD.vcf.gz; # fi;
tabix ${vcf/.PL.vcf.gz/.PL.AD.vcf.gz} -f
