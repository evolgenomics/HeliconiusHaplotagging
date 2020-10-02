#!/bin/bash

#$ -l h_vmem=8G
#$ -t 57-196
#$ -N mpileup
#$ -o /fml/chones/projects/PC062_Heliconius_Haplo/STITCH/job_out/
#$ -j y
#$ -S /bin/bash
#$ -cwd

if [ ! -e "/tmp/frankyc/PC062_Heliconius_Haplo/mpileup" ];then mkdir -p /tmp/frankyc/PC062_Heliconius_Haplo/mpileup;fi
#
#intervals=$1;
#task=`printf %04d $SGE_TASK_ID`;
##task=`printf %04d $2`;
#echo $intervals" - "$task
#
#fbname=$(basename $intervals .intervals)
#echo $fbname
#
dir=`pwd`
#echo $dir
#
region=`awk 'NR == "'$SGE_TASK_ID'"' helera1_demo.chr.intervals`;
output=PC062_helera1_demo_merged.mpileup.$region.vcf
chr=`echo $region | sed 's/:1-.\+//'`
pos=STITCH/helera1_demo_pos/PC062_helera1_demo_merged.Groups.$chr.QUAL20_SNP_REF_2_ALT_2.pos

/fml/chones/local/bin/samtools mpileup -r $region -a -l $pos -g -t AD,DP,DV,DPR,INFO/DPR,DP4,SP -q 10 -Q 20 -R `cut -c 4- STITCH/bams_samples_helera1_demo.list` -f /fml/chones/genome/gbdb/helera1_demo/hel
era1_demo.fa -r $region -F 1024 > /tmp/frankyc/PC062_Heliconius_Haplo/mpileup/$output

mv /tmp/frankyc/PC062_Heliconius_Haplo/mpileup/$output $dir/STITCH/out_helera1_demo/
