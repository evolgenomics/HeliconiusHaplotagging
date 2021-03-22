#!/bin/bash

#Example: ../hapcutVcf.helera1_demo.sh hap_0.linkedReads.diploid.noIndel.bam stitch.Herato1603.450000.550000.PL.vcf.gz noIndel
#Example2:  ../hapcutVcf.helera1_demo.sh hap_0.linkedReads.diploid.linked.noIndel.bam stitch.Herato1603.450000.550000.PL.vcf.gz linked.noIndel

bam=$1
vcf=$2
tag=$3
runBC=`basename $bam | sed 's/\./\t/g' | cut -f 1`;
runBC_base=`basename $bam .bam`
sampleHetVCF=$runBC.$tag.PL.het.vcf
sampleVCF=$runBC.$tag.PL.vcf

out=`basename $bam .$tag.bam`
out_dir=tmp
tmp=/tmp/$runBC_base
if [ ! -e "$tmp" ]; then mkdir -p $tmp; fi

echo "BAM: $bam"
echo "vcf: $vcf"
echo "runBC: $runBC"
echo "tag: $tag"
echo "sampleHetVCF: $sampleHetVCF"
echo "sampleVCF: $sampleVCF"
rm $tmp/$out.$tag.*

echo "Writing temporary VCF files..." 
bcftools view -s $runBC $vcf -i 'INFO/INFO_SCORE >= 0.2' | awk '/^#/;/CHROM/ {OFS="\t"}; !/^#/ && $10~/^0\/1/' > $tmp/$sampleHetVCF
#bcftools view -s $runBC $vcf -i 'INFO/INFO_SCORE >= 0.2'  | awk '/^#/;/CHROM/ {OFS="\t"}; !/^#/ &&  $10~/^0\/0/ {$10="0|0:"substr($10,5);print $0};  !/^#/ && $10~/^0\/1/; !/^#/ &&  $10~/^1\/1/ {$10="1|1:"substr($10,5);print $0}; !/^#/ {print $0}' > $tmp/$sampleVCF
bcftools view -s $runBC $vcf -i 'INFO/INFO_SCORE >= 0.2'  | awk '/^#/;/CHROM/ {OFS="\t"}; !/^#/ &&  $10~/^0\/0/ {$10="0|0:"substr($10,5);print $0}; !/^#/ &&  $10~/^1\/1/ {$10="1|1:"substr($10,5);print $0}; !/^#/ && $10~/^0\/1/ {$9=substr($9, 4); $10=substr($10,5);print $0}' > $tmp/$sampleVCF
#Option to use if PG field has not been used as backup GT
#bcftools view -s $runBC $vcf -i 'INFO/INFO_SCORE >= 0.2' | awk '/^#/;/CHROM/ {OFS="\t"}; !/^#/ && $10~/^0\/1/; !/^#/ && $10~/\.\/\./ {split($10,FMT,":"); split(FMT[2],GP,",");if (GP[2] > GP[1] && GP[2] > GP[3]) {$10="0/1"substr($10,4);print $0}}' > $tmp/$sampleHetVCF
#bcftools view -s $runBC $vcf -i 'INFO/INFO_SCORE >= 0.2'  | awk '/^#/;/CHROM/ {OFS="\t"}; !/^#/ &&  $10~/^0\/0/ {$10=$10":0|0";$9=$9":PG";print $0};  !/^#/ && $10~/^0\/1/; !/^#/ &&  $10~/^1\/1/ {$10=$10":1|1";$9=$9":PG";print $0}; !/^#/ && $10~/\.\/\./ {split($10,FMT,":"); split(FMT[2],GP,",");if (GP[1] > GP[2] && GP[1] > GP[3]) {$10=$10":0|0";$9=$9":PG";print $0}; if(GP[3]>GP[2] && GP[3] > GP[1]) {$10=$10":1|1";$9=$9":PG";print $0};if (GP[2] > GP[1] && GP[2] > GP[3]) {print $0}; if (GP[1] == GP[2] || GP[1] == GP[3] || GP[2] == GP[3]){print $0} }' > $tmp/$sampleVCF
# 

##
echo "Running HAPCUT2..." 
samtools index ../$bam;
extractHAIRS --bam ../$bam --VCF $tmp/$sampleHetVCF --out $tmp/$out.$tag.raw.fragments 
extractHAIRS --10X 1 --bam ../$bam --VCF $tmp/$sampleHetVCF --out $tmp/$out.$tag.unlinked.fragments 
grep -ax '.*' $tmp/$out.$tag.unlinked.fragments > $tmp/$out.$tag.unlinked.fragments.1; 
mv $tmp/$out.$tag.unlinked.fragments.1 $tmp/$out.$tag.unlinked.fragments;
python3 ${HAPCUT2}/utilities/LinkFragments.py  --bam ../$bam --VCF $tmp/$sampleHetVCF --fragments $tmp/$out.$tag.unlinked.fragments --out $tmp/$out.$tag.linked.fragments -d 50000;
##
HAPCUT2 --fragments $tmp/$out.$tag.linked.fragments --VCF $tmp/$sampleHetVCF --out $tmp/$out.$tag.hapcut2.threshold_30.output --nf 1 --threshold 30 --error_analysis_mode 1 --call_homozygous 1 --outvcf 1   
HAPCUT2 --fragments $tmp/$out.$tag.raw.fragments --VCF $tmp/$sampleHetVCF --out $tmp/$out.$tag.hapcut2.threshold_30.unlinked.output --threshold 30 --error_analysis_mode 1 --call_homozygous 1 --outvcf 1   
###
echo "Processing HAPCUT2 output..." 
../../scripts/hapCUT2_to_bed.sh $tmp/$out.$tag.hapcut2.threshold_30.output
../../scripts/hapCUT2_to_bed.sh $tmp/$out.$tag.hapcut2.threshold_30.unlinked.output
###
../../scripts/n50_extract.sh $tmp/$out.$tag.hapcut2.threshold_30.bed
../../scripts/n50_extract.sh $tmp/$out.$tag.hapcut2.threshold_30.unlinked.bed
###
echo "Creating annotation files..." 
bcftools query -f "%CHROM\t%POS[\t%GT\t%PS\t%PQ\t%PD]\n" $tmp/$out.$tag.hapcut2.threshold_30.output.phased.VCF | bgzip -c > $tmp/$out.$tag.hapcut2.threshold_30.output.phased.annot.gz
bcftools query -f "%CHROM\t%POS[\t%GT\t%PS\t%PQ\t%PD]\n" $tmp/$out.$tag.hapcut2.threshold_30.unlinked.output.phased.VCF | bgzip -c > $tmp/$out.$tag.hapcut2.threshold_30.unlinked.output.phased.annot.gz
###
tabix -b 2 -e 2 $tmp/$out.$tag.hapcut2.threshold_30.output.phased.annot.gz
tabix -b 2 -e 2 $tmp/$out.$tag.hapcut2.threshold_30.unlinked.output.phased.annot.gz
###
#echo "bcftools annotate -h ../../add.hdr -a $out.$tag.hapcut2.threshold_30.output.phased.annot.gz $sampleVCF -c CHROM,POS,FMT/PG,FMT/PS,FMT/PQ,FMT/PD -m +HAPCUT=1 | bcftools view - -Oz -o ${sampleVCF/.vcf/.HAPCUT2.vcf.gz};" 
##echo "Combining annotations into single VCF..." 
###bcftools annotate -h add.hdr -a $tmp/$out.$tag.hapcut2.threshold_30.output.phased.annot.gz $tmp/$sampleVCF -c CHROM,POS,FMT/GT,FMT/PS,FMT/PQ,FMT/PD -m +HAPCUT=1 | /fml/chones/local/bin/bcftools view - -Oz -o $tmp/${sampleVCF/.vcf/.HAPCUT2.vcf.gz}; 
bcftools annotate -h ../add.hdr -a $tmp/$out.$tag.hapcut2.threshold_30.output.phased.annot.gz $tmp/$sampleVCF -c CHROM,POS,FMT/GX,FMT/PS,FMT/PQ,FMT/PD -m +HAPCUT=1 |  awk '!/<ID=GX/' | sed 's/:GX:/:GT:/' | /fml/chones/local/bin/bcftools view - -Oz -o $tmp/${sampleVCF/.vcf/.HAPCUT2.vcf.gz}; 
bcftools annotate -h ../add.hdr -a $tmp/$out.$tag.hapcut2.threshold_30.unlinked.output.phased.annot.gz $tmp/$sampleVCF -c CHROM,POS,FMT/GX,FMT/PS,FMT/PQ,FMT/PD -m +HAPCUT=1 |  awk '!/<ID=GX/' | sed 's/:GX:/:GT:/' | /fml/chones/local/bin/bcftools view - -Oz -o $tmp/${sampleVCF/.vcf/.HAPCUT2_unlinked.vcf.gz}; 
tabix $tmp/${sampleVCF/.vcf/.HAPCUT2.vcf.gz}
tabix $tmp/${sampleVCF/.vcf/.HAPCUT2_unlinked.vcf.gz}
##
##rm $tmp/$sampleHetVCF
##rm $tmp/$sampleVCF
#rm $tmp/$out.$tag.unlinked.fragments
#rm $tmp/$out.$tag.linked.fragments
#rm $tmp/$out.$tag.hapcut2.threshold_30.*output.phased.annot.gz*

mv -v $tmp/${sampleVCF/.vcf/.HAPCUT2.vcf.gz}* $out_dir/
mv -v $tmp/${sampleVCF/.vcf/.HAPCUT2_unlinked.vcf.gz}* $out_dir/
mv -v $tmp/$out.$tag.hapcut2.threshold_30*output $out_dir/
mv -v $tmp/$out.$tag.hapcut2.threshold_30.*bed $out_dir/
