#!/bin/bash

#$ -l h_vmem=2G
#$ -t 1-264
#$ -N struct_jac_redo
#$ -o /fml/chones/projects/PC062_Heliconius_Haplo/struct/job_out/
#$ -j y
#$ -S /bin/bash
#$ -cwd

if [ ! -e "/tmp/frankyc/PC062_Heliconius_Haplo/struct" ];then mkdir -p /tmp/frankyc/PC062_Heliconius_Haplo/struct;fi

intervals=$1;
task=`printf %04d $SGE_TASK_ID`;
#task=`printf %04d $2`;
echo $intervals" - "$task

fbname=$(basename $intervals .intervals)
echo $fbname

dir=$(dirname $intervals)
echo $dir


perl SCRIPTS/structural_jaccard.hmel2.5.1.pl joblist/$intervals.chunk_$task.list /tmp/frankyc/PC062_Heliconius_Haplo/struct/$fbname.interval_$task.redo_sites.jaccard_matrix 

mv /tmp/frankyc/PC062_Heliconius_Haplo/struct/$fbname.interval_$task.redo_sites.jaccard_matrix  $dir/out/
rm /tmp/frankyc/PC062_Heliconius_Haplo/struct/$fbname.interval_$task.jaccard_matrix 
