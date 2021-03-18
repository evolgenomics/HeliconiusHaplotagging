#!/bin/bash
echo "
[setup.sh] Creating folder structure...

";
for folder in sourceData altFa simARG simTrees simBam simFastq simHaps simMuts log; do \
	mkdir $folder;
done

chmod 755 scripts/*

echo "
[setup.sh] Getting and indexing the reference genome assembly...

";
cd sourceData
wget http://ftp.tuebingen.mpg.de/fml/ag-chan/haplotagging/simulations/helera1_demo_Herato1603.fa 2&> ../log/setup.err.log
bwa index helera1_demo_Herato1603.fa
samtools faidx helera1_demo_Herato1603.fa

echo "
[setup.sh] Downloading the reference phased VCF and generating the source haplotypes...

";
wget http://ftp.tuebingen.mpg.de/fml/ag-chan/haplotagging/simulations/PC062_merged_Herato1603.3.45Mb.PL.AD.HAPCUT2.vcf.gz 2&> ../log/setup.err.log 
wget http://ftp.tuebingen.mpg.de/fml/ag-chan/haplotagging/simulations/PC062_merged_Herato1603.3.45Mb.PL.AD.HAPCUT2.vcf.gz.tbi 2&> ../log/setup.err.log 
wget http://ftp.tuebingen.mpg.de/fml/ag-chan/haplotagging/simulations/add.hdr 2&> ../log/setup.err.log 
wget http://ftp.tuebingen.mpg.de/fml/ag-chan/haplotagging/simulations/header_cols.vcf 2&> ../log/setup.err.log 
wget http://ftp.tuebingen.mpg.de/fml/ag-chan/haplotagging/simulations/vcf.header 2&> ../log/setup.err.log 

bcftools query -f "[\t%GT]\n" PC062_merged_Herato1603.3.45Mb.PL.AD.HAPCUT2.vcf.gz | sed 's/[|/]/\t/g;s/^\t//' | datamash transpose | sed 's/\t//g' | sort > PC062_merged_Herat1603_3.45Mb.simBlock.horiz.sorted 
bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\n" PC062_merged_Herato1603.3.45Mb.PL.AD.HAPCUT2.vcf.gz > PC062_merged_Herato1603.3.45Mb.PL.AD.HAPCUT2.pos
for i in `seq 0 2 967`; do echo "hap_"$i; done > sample.names
echo "Herato1603	3350000	3650000" > focal_region.bed 
cd ..
