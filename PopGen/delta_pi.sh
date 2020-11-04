#!/bin/bash

module load htslib/1.2.1 bcftools-1.9-gcc-5.4.0-b2hdt5n samtools/1.3.1 vcftools

# Compute 10 kb pi for each scaffold
for file in ../vcf/run164_merged_*.PL.AD.HAPCUT2.inf0.05.vcf.gz
do 
  prefix=`basename $i`
  prefix=${prefix%.vcf.gz}
  vcftools --gzvcf $file --keep malleti.samples --window-pi 10000 --window-pi-step 2500 --out $prefix.malleti.10kb
  vcftools --gzvcf $file --keep plesseni.samples --window-pi 10000 --window-pi-step 2500 --out $prefix.plesseni.10kb
done

# Compute 50 kb pi for each scaffold
for file in ../vcf/run164_merged_*.PL.AD.HAPCUT2.inf0.05.vcf.gz
do 
  prefix=`basename $i`
  prefix=${prefix%.vcf.gz}
  vcftools --gzvcf $file --keep malleti.samples --window-pi 50000 --out $prefix.malleti.50kb
  vcftools --gzvcf $file --keep plesseni.samples --window-pi 50000 --out $prefix.plesseni.50kb
done

# Delta pi was then computed in R as pi(highland subspecies) - pi(lowland subspecies)
