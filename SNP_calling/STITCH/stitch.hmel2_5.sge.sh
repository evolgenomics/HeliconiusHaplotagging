#!/bin/bash

#$ -l h_vmem=16G
#$ -N stitch_hmel2.5
#$ -o /fml/chones/projects/PC062_Heliconius_Haplo/STITCH/job_out/
#$ -j y
#$ -S /bin/bash
#$ -cwd

if [ ! -e "/tmp/frankyc/PC062_Heliconius_Haplo/stitch" ];then mkdir -p /tmp/frankyc/PC062_Heliconius_Haplo/stitch;fi

chr=`sed 's/:/\t/;s/-/\t/' SCRIPT/hmel2.5.intervals | awk 'NR == '$SGE_TASK_ID' {print $1}'` 
start=`sed 's/:/\t/;s/-/\t/' SCRIPT/hmel2.5.intervals | awk 'NR == '$SGE_TASK_ID' {print $2+1}'` 
end=`sed 's/:/\t/;s/-/\t/' SCRIPT/hmel2.5.intervals | awk 'NR == '$SGE_TASK_ID' {print $3}'` 
pos=hmel2.5_pos/Run164_merged.$chr.QUAL20_SNP_REF_2_ALT_2.refAlt.pos 
echo $SGE_TASK_ID" "$chr" "$start" "$end
/fml/chones/local/STITCH/STITCH.R --chr=$chr --regionStart=$start --regionEnd=$end --buffer=25000 --bamlist=bams_samples_hmel.list --sampleNames_file=bams_samples_hmel.names.list --posfile=$pos --outputdir=/tmp/frankyc/PC062_Heliconius_Haplo/stitch/stitch_hmel2.5_${chr}_${start}_$end --K=40 --nCores=1 --method=diploid --nGen=500 --readAware=TRUE --keepInterimFiles=FALSE --iSizeUpperLimit=500000 --keepSampleReadsInRAM=TRUE
mv -v /tmp/frankyc/PC062_Heliconius_Haplo/stitch/stitch_hmel2.5_${chr}_${start}_${end}* out/
