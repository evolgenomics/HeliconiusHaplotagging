#!/bin/bash

#$ -l h_vmem=8G
#$ -t 1-76
#$ -N omegaPlus
#$ -o /fml/chones/projects/PC062_Heliconius_Haplo/split_vcfs/job_out/
#$ -j y
#$ -S /bin/bash
#$ -cwd

if [ ! -e "/tmp/frankyc/PC062_Heliconius_Haplo/omega" ];then mkdir -p /tmp/frankyc/PC062_Heliconius_Haplo/omega;fi

species=`awk 'NR =="'$SGE_TASK_ID'" {print $1}' SCRIPTS/omega.melpomene.list`
vcf=`awk 'NR =="'$SGE_TASK_ID'" {print $2}' SCRIPTS/omega.melpomene.list`

echo $species"	"$vcf
cat MaCS.in.template <(/fml/chones/local/bin/bcftools query -f "%POS\t[%GT]\n" -S $species.samples $vcf | sed 's/.\/./0|1/g;s/|//g'| awk 'NR ==1 {s = length($2)};{print "SITE:\t"NR-1"\t0."(sprintf("%08d",$1))"\t"$2}; END {print "TOTAL_SAMPLES:    "s;print "TOTAL_SITES:\t"NR;print "BEGIN_SELECTED_SITES";str="";for(i=0;i<NR;i++){str=str"\t"i};print NR >"/tmp/frankyc/PC062_Heliconius_Haplo/omega/'$SGE_TASK_ID'.tmp"; print substr(str,2)}') > /tmp/frankyc/PC062_Heliconius_Haplo/omega/${vcf/.vcf.gz/.focal_$species.MaCS.in};
cd /tmp/frankyc/PC062_Heliconius_Haplo/omega/;
/fml/chones/local/bin/OmegaPlus -name ${vcf/.vcf.gz/.focal_$species} -input /tmp/frankyc/PC062_Heliconius_Haplo/omega/${vcf/.vcf.gz/.focal_$species.MaCS.in} -grid `tail -5 /tmp/frankyc/PC062_Heliconius_Haplo/omega/${vcf/.vcf.gz/.focal_$species.MaCS.in} | head -1 | awk '{print sprintf("%.0f", ($3 * 100000000/10))}'` -maxwin 500 -minwin 5 -binary -length 100000000 -seed 3141323 -all -impute N
mv -v /tmp/frankyc/PC062_Heliconius_Haplo/omega/Omega*.${vcf/.vcf.gz/.focal_$species} /fml/chones/projects/PC062_Heliconius_Haplo/split_vcfs/omega 
cd /fml/chones/projects/PC062_Heliconius_Haplo/split_vcfs/
