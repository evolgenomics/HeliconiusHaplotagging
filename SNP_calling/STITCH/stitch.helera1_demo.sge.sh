#!/bin/bash

#$ -l h_vmem=16G
#$ -t 1-943
#$ -N stitch_helera1
#$ -o /fml/chones/projects/PC062_Heliconius_Haplo/STITCH/job_out/
#$ -j y
#$ -S /bin/bash
#$ -cwd

if [ ! -e "/tmp/frankyc/PC062_Heliconius_Haplo/stitch" ];then mkdir -p /tmp/frankyc/PC062_Heliconius_Haplo/stitch;fi

chr=`sed 's/:/\t/;s/-/\t/' SCRIPT/helera1_demo.intervals.1 | awk 'NR == '$SGE_TASK_ID' {print $1}'` 
start=`sed 's/:/\t/;s/-/\t/' SCRIPT/helera1_demo.intervals.1 | awk 'NR == '$SGE_TASK_ID' {print $2+1}'` 
end=`sed 's/:/\t/;s/-/\t/' SCRIPT/helera1_demo.intervals.1 | awk 'NR == '$SGE_TASK_ID' {print $3}'` 
pos=helera1_demo_pos/PC062_helera1_demo_merged.Groups.$chr.QUAL20_SNP_REF_2_ALT_2.pos 
echo $SGE_TASK_ID" "$chr" "$start" "$end
/fml/chones/local/STITCH/STITCH.R --chr=$chr --regionStart=$start --regionEnd=$end --buffer=25000 --bamlist=bams_samples_helera1_demo.list --sampleNames_file=bams_samples_helera1_demo.names.list --posfile=$pos --outputdir=/tmp/frankyc/PC062_Heliconius_Haplo/stitch/stitch_helera1_demo_${chr}_${start}_$end --K=30 --method=diploid --nGen=500 --nCores=1 --readAware=TRUE --keepInterimFiles=FALSE --shuffle_bin_radius=500 --expRate=5 --iSizeUpperLimit=500000 --keepSampleReadsInRAM=TRUE

mv -v /tmp/frankyc/PC062_Heliconius_Haplo/stitch/stitch_helera1_demo_${chr}_${start}_$end out_helera1_demo/stitch_helera1_demo_${chr}_${start}_$end
