#!/bin/bash

tag=$1;

mkdir ${tag}_subsample
#cd ${tag}_subsample;
#for i in `seq 0 2 967`; do samtools view hap_$i.linkedReads.diploid.$tag.bam | sed 's/BX:Z:/\nBX:Z:/' | grep BX:Z | sort | uniq > ${tag}_subsample/hap_$i.linkedReads.diploid.$tag.bx; done
for frac in 0.01 0.025 0.05 0.075 0.1 0.25 0.5 0.75; do 
#for frac in 0.1 0.25 0.5 0.75; do 
	echo "Subsampling at ${frac}x...";
	perl mol_subsample.pl $tag $frac
	#for i in `seq 0 2 967`; do 
		
		#molNum=`cat hap_$i.linkedReads.diploid.$tag.bx | wc -l`; 
		#num=`echo $molNum | awk '{print sprintf("%.0f", $1*'$frac')}'`;
		#echo "^@" > hap_$i.linkedReads.diploid.$tag.$frac.grep;
		#shuf -n $num hap_$i.linkedReads.diploid.$tag.bx | awk '{print $0"\\b"}' >> hap_$i.linkedReads.diploid.$tag.$frac.grep;
		#samtools view ../hap_$i.linkedReads.diploid.$tag.bam -h | grep -f hap_$i.linkedReads.diploid.$tag.$frac.grep | samtools view - -O BAM -o hap_$i.linkedReads.diploid.$tag.$frac.bam;
		#samtools view ../hap_$i.linkedReads.diploid.linked.$tag.bam -h | grep -f hap_$i.linkedReads.diploid.$tag.$frac.grep | samtools view - -O BAM -o hap_$i.linkedReads.diploid.linked.$tag.$frac.bam;
		#samtools index -@ 32 hap_$i.linkedReads.diploid.$tag.$frac.bam;
		#samtools index -@ 32 hap_$i.linkedReads.diploid.linked.$tag.$frac.bam;
	#done;
	ls ${tag}_subsample/hap_$i.linkedReads.diploid.$tag.$frac.bam | sort -k 1,1V > bamfiles.linkedReads.$tag.$frac.968.list
	ls ${tag}_subsample/hap_$i.linkedReads.diploid.linked.$tag.$frac.bam | sort -k 1,1V > bamfiles.linkedReads.hacked.$tag.$frac.968.list
	#rm ${tag}_subsample/hap_*.linkedReads.diploid.$tag.$frac.grep
done
rm ${tag}_subsample/hap*.linkedReads.diploid.$tag.bx
cd ..
