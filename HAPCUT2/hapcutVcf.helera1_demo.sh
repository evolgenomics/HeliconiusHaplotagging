#!/bin/bash

#This script is used to create a combined, phased VCF file for a single individual sample
#It does so with the following steps:
#
#1. generate two temporary VCF files - one with only heterozygous sites, and the other one with all sites
#2. HAPCUT2 pipeline, with basically three steps:
#2a. Parse the bam file for reads that correspond to the heterozygous sites for HAPCUT2 using the extractHAIRs utility with the 10X flag on to generate a BX-tagged, unlinked fragment file
#2b. Link the fragments using HAPCUT2's LinkFragments.py utility.
#2c. Run HAPCUT2 proper, with VCF output option
#3. Perform some basic data extraction from the HAPCUT2 output
#4. Merge the full VCF file with the HAPCUT2-derived VCF files and annotating the INFO and FORMAT fields accordingly
#5. Clean-up

bam=$1
vcf=$2
region=$3
runBC=`basename $bam .bam`
camid=`echo $runBC | sed 's/Run/run/;s/demo.split.//;s/.merged//;s/.bam$//' | xargs -i grep {} ../SAMPLE_INFO.barcodes.list | cut -f 4 | head -1`
sampleHetVCF=$runBC.pMarkdup.$region.PL.AD.het.vcf
sampleVCF=$runBC.pMarkdup.$region.PL.AD.vcf
chrom=`echo $region | sed 's/:/\t/' | cut -f 1`;

out=`basename $bam .bam`
out_dir=/fml/chones/projects/PC062_Heliconius_Haplo/split_vcfs/$chrom
tmp=/tmp/frankyc/PC062_Heliconius_Haplo/fragments/$camid
if [ ! -e "$tmp" ]; then mkdir -p $tmp; fi

echo "BAM: $bam"
echo "vcf: $vcf"
echo "runBC: $runBC"
echo "camID: $camid"
echo "chrom: $chrom"
echo "region: $region"
echo "sampleHetVCF: $sampleHetVCF"
echo "sampleVCF: $sampleVCF"
rm $tmp/$out.$region.*

echo "Writing temporary VCF files..." 
/fml/chones/local/bin/bcftools view -s $camid $vcf -i 'INFO/INFO_SCORE >= 0.2' | awk '/^#/;/CHROM/ {OFS="\t"}; !/^#/ && $10~/^0\/1/' > $tmp/$sampleHetVCF
#/fml/chones/local/bin/bcftools view -s $camid $vcf -i 'INFO/INFO_SCORE >= 0.2'  | awk '/^#/;/CHROM/ {OFS="\t"}; !/^#/ &&  $10~/^0\/0/ {$10="0|0:"substr($10,5);print $0};  !/^#/ && $10~/^0\/1/; !/^#/ &&  $10~/^1\/1/ {$10="1|1:"substr($10,5);print $0}; !/^#/ {print $0}' > $tmp/$sampleVCF
/fml/chones/local/bin/bcftools view -s $camid $vcf -i 'INFO/INFO_SCORE >= 0.2'  | awk '/^#/;/CHROM/ {OFS="\t"}; !/^#/ &&  $10~/^0\/0/ {$10="0|0:"substr($10,5);print $0}; !/^#/ &&  $10~/^1\/1/ {$10="1|1:"substr($10,5);print $0}; !/^#/ && $10~/^0\/1/ {$9=substr($9, 4); $10=substr($10,5);print $0}' > $tmp/$sampleVCF
#Option to use if PG field has not been used as backup GT
#/fml/chones/local/bin/bcftools view -s $camid $vcf -i 'INFO/INFO_SCORE >= 0.2' | awk '/^#/;/CHROM/ {OFS="\t"}; !/^#/ && $10~/^0\/1/; !/^#/ && $10~/\.\/\./ {split($10,FMT,":"); split(FMT[2],GP,",");if (GP[2] > GP[1] && GP[2] > GP[3]) {$10="0/1"substr($10,4);print $0}}' > $tmp/$sampleHetVCF
#/fml/chones/local/bin/bcftools view -s $camid $vcf -i 'INFO/INFO_SCORE >= 0.2'  | awk '/^#/;/CHROM/ {OFS="\t"}; !/^#/ &&  $10~/^0\/0/ {$10=$10":0|0";$9=$9":PG";print $0};  !/^#/ && $10~/^0\/1/; !/^#/ &&  $10~/^1\/1/ {$10=$10":1|1";$9=$9":PG";print $0}; !/^#/ && $10~/\.\/\./ {split($10,FMT,":"); split(FMT[2],GP,",");if (GP[1] > GP[2] && GP[1] > GP[3]) {$10=$10":0|0";$9=$9":PG";print $0}; if(GP[3]>GP[2] && GP[3] > GP[1]) {$10=$10":1|1";$9=$9":PG";print $0};if (GP[2] > GP[1] && GP[2] > GP[3]) {print $0}; if (GP[1] == GP[2] || GP[1] == GP[3] || GP[2] == GP[3]){print $0} }' > $tmp/$sampleVCF
# 
export LD_LIBRARY_PATH=/fml/chones/local/lib
export PYTHONPATH=/fml/chones/local/lib/python3.6/site-packages:$PYTHONPATH
#
echo "Running HAPCUT2..." 
/fml/chones/local/bin/extractHAIRS --10X 1 --bam $bam --VCF $tmp/$sampleHetVCF --region $region --out $tmp/$out.$region.unlinked.fragments 
grep -ax '.*' $tmp/$out.$region.unlinked.fragments > $tmp/$out.$region.unlinked.fragments.1; 
mv $tmp/$out.$region.unlinked.fragments.1 $tmp/$out.$region.unlinked.fragments;
python3 /fml/chones/local/HapCUT2/utilities/LinkFragments.py  --bam $bam --VCF $tmp/$sampleHetVCF --fragments $tmp/$out.$region.unlinked.fragments --out $tmp/$out.$region.linked.fragments -d 50000;
#
/fml/chones/local/bin/HAPCUT2 --fragments $tmp/$out.$region.linked.fragments --VCF $tmp/$sampleHetVCF --out $tmp/$out.$region.hapcut2.threshold_30.output --nf 1 --threshold 30 --error_analysis_mode 1 --call_homozygous 1 --outvcf 1   
##
#echo "Processing HAPCUT2 output..." 
#/fml/chones/projects/PC062_Heliconius_Haplo/hapcut2/hapCUT2_to_bed.sh $tmp/$out.$region.hapcut2.threshold_30.output
##
#/fml/chones/projects/PC062_Heliconius_Haplo/hapcut2/n50_extract.sh $tmp/$out.$region.hapcut2.threshold_30.bed
##
#echo "Creating annotation files..." 
/fml/chones/local/bin/bcftools query -f "%CHROM\t%POS[\t%GT\t%PS\t%PQ\t%PD]\n" $tmp/$out.$region.hapcut2.threshold_30.output.phased.VCF | bgzip -c > $tmp/$out.$region.hapcut2.threshold_30.output.phased.annot.gz
##
tabix -b 2 -e 2 $tmp/$out.$region.hapcut2.threshold_30.output.phased.annot.gz
##
echo "bcftools annotate -h add.hdr -a $out.$region.hapcut2.threshold_30.output.phased.annot.gz $sampleVCF -c CHROM,POS,FMT/PG,FMT/PS,FMT/PQ,FMT/PD -m +HAPCUT=1 | bcftools view - -Oz -o ${sampleVCF/.vcf/.HAPCUT2.vcf.gz};" 
echo "Combining annotations into single VCF..." 
##/fml/chones/local/bin/bcftools annotate -h add.hdr -a $tmp/$out.$region.hapcut2.threshold_30.output.phased.annot.gz $tmp/$sampleVCF -c CHROM,POS,FMT/GT,FMT/PS,FMT/PQ,FMT/PD -m +HAPCUT=1 | /fml/chones/local/bin/bcftools view - -Oz -o $tmp/${sampleVCF/.vcf/.HAPCUT2.vcf.gz}; 
/fml/chones/local/bin/bcftools annotate -h add.hdr -a $tmp/$out.$region.hapcut2.threshold_30.output.phased.annot.gz $tmp/$sampleVCF -c CHROM,POS,FMT/GX,FMT/PS,FMT/PQ,FMT/PD -m +HAPCUT=1 |  awk '!/<ID=GX/' | sed 's/:GX:/:GT:/' | /fml/chones/local/bin/bcftools view - -Oz -o $tmp/${sampleVCF/.vcf/.HAPCUT2.vcf.gz}; 
tabix $tmp/${sampleVCF/.vcf/.HAPCUT2.vcf.gz}
#
echo "Cleaning up"
rm $tmp/$sampleHetVCF
rm $tmp/$sampleVCF
rm $tmp/$out.$region.unlinked.fragments
rm $tmp/$out.$region.linked.fragments
rm $tmp/$out.$region.hapcut2.threshold_30.output.phased.annot.gz*

mv -v $tmp/${sampleVCF/.vcf/.HAPCUT2.vcf.gz}* $out_dir/
mv -v $tmp/$out.$region.hapcut2.threshold_30.output $out_dir/
mv -v $tmp/$out.$region.hapcut2.threshold_30.*bed $out_dir/
