#!/bin/bash

#$ -pe parallel 1
#$ -l h_vmem=15G
#$ -N MarkDupl
#$ -o /fml/mickle/data/Marek/Jobs/
#$ -j y
#$ -S /bin/bash
#$ -cwd

if [ ! -e "/tmp/mkucka" ];then mkdir /tmp/mkucka;fi


bam=$1
vcf=$2
chr=$3
out=`basename $bam .bam`

#unset PYTHONPATH
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}
export PYTHONPATH=${PYTHONPATH}

extractHAIRS --10X 1 --bam $bam --VCF $vcf --out $out.$chr.unlinked.fragments
awk '!/[ABCD]00/' $out.$chr.unlinked.fragments > $out.$chr.no00.unlinked.fragments

python3 /fml/chones/local/HapCUT2/utilities/LinkFragments.py  --bam $bam --VCF $vcf --fragments $out.$chr.no00.unlinked.fragments --out $out.$chr.no00.linked.fragments -d 50000;

HAPCUT2 --fragments $out.$chr.no00.linked.fragments --VCF $vcf --out $out.Chr$chr.hapcut2.threshold_30.errAn.output --nf 1 --threshold 30
