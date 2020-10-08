mpileup.sh
#!/bin/bash
input=$1
bed=$2
samtools mpileup -g -t INFO/AD,AD,DP,DV,DPR,INFO/DPR,DP4,SP -q 10 -Q 20 -f /fml/chones/genome/gbdb/mm10/mm10.fa -R -l $bed -F 1024 $input | bcftools call -f GQ,GP -O v -v -m -o $bed.vcf 
