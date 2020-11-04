#!/bin/bash

# Script to run NGSadmix with all H. melpomene samples

module load htslib/1.2.1 bcftools-1.9-gcc-5.4.0-b2hdt5n samtools/1.3.1 vcftools

# Filter the vcf files a bit more stringently (INFO score at least 0.5 instead of 0.2)
file=run164_merged_${scaff}.PL.AD.HAPCUT2.vcf.gz
bcftools view -i "INFO_SCORE >= 0.5"  -Oz -o ${file%.vcf.gz}.inf0.05.vcf.gz $file

# Randomly subsample sites to 10% to decrease computing time and remove LD 
file=${file%.vcf.gz}.inf0.05
zcat $file | sed 's/ID=HAPCUT=1/ID=HAPCUT/' | vcfrandomsample -r 0.1 | gzip > ${file%.vcf.gz}.subsampled.vcf.gz

# Concatenate the vcf files of all scaffolds to a single file
bcftools concat -Oz -o run164_merged_PL.AD.HAPCUT2.inf0.05.subsampled.vcf.gz run164_merged_Hmel2*.PL.AD.HAPCUT2.inf0.05.subsampled.vcf.gz

# Fix missing genotypes:
file=run164_merged_PL.AD.HAPCUT2.inf0.05.subsampled
zcat $file.vcf.gz | sed -e 's|./.:.:.:.:.:.|0/0:./.:0.33,0.33,0.34:0:10,10,10:0,0|g' | gzip > $file.corr.vcf.gz

# Generate a Beagle file with angsd
file=run164_merged_PL.AD.HAPCUT2.inf0.05.subsampled.corr
angsd -out $file -nThreads 6 \
  -vcf-gl $file.vcf.gz -doGlf 2 \
  -doMajorMinor 1 -doMaf 2 -SNP_pval 1e-6 

# Run NGSadmix with different K values
# run it 10 times per K value to assess the best number of clusters (K) with Clumpak
file=run164_merged_PL.AD.HAPCUT2.inf0.05.subsampled.corr
for j in {1..10}
do
    # K ranging from 1 to 4
    for k in {1..4}
    do
        NGSadmix -likes ${file}.beagle.gz -minMaf 0.05 -K $k -o ${file}.K${k}.run${j} 
    done
done

# Collect the likelihoods:
for log in ls *.log; do echo $log >> logfile ; grep -Po 'like=\K[^ ]+' $log >> logfile; done

# Reformat for Clumpak:
cat logfile | tr '\n' ' ' | sed -e 's/ stitch_hmel2.5.PL.AD.inf0.5.subsampled./\n/g' -e 's/.log//g' -e 's/\.run.[0-9]*//g' \
 -e 's/stitch_hmel2.5.PL.AD.inf0.5.subsampled.//g' -e 's/K//g' -e 's/ /\t/g' | grep -v maf > likelihoods.forclumpak

# Load the likelihoods.forclumpak file as “Log probability table” to http://clumpak.tau.ac.il/bestK.html. (http://clumpak.tau.ac.il/bestK.html)

# Get 2 columns of the K=3 run for GWAS 
cut -f 1,2 -d " " $file.maf0.05.K3.qopt > NGSadmix.K3.prop
