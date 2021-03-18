#!/bin/bash
tag=$1;
folder=$2;

cd $folder;
#../add_PL.helera1_demo.sh stitch.Herato1603.450000.550000.vcf.gz
#bcftools reheader -s ../sample.names stitch.Herato1603.450000.550000.PL.vcf.gz > stitch.Herato1603.450000.550000.PL.reheadered.vcf.gz
#bcftools index stitch.Herato1603.450000.550000.PL.reheadered.vcf.gz;
#perl ../check_truth_ARG.pl $tag
#mkdir tmp;


for i in `seq 0 2 60`; do echo "HAPCUT - $i"; ../hapcutVcf.helera1_demo.sh hap_$i.linkedReads.diploid.$tag.bam stitch.Herato1603.450000.550000.PL.reheadered.vcf.gz $tag; done
bcftools merge -Oz -o stitch.Herato1603.450000.550000.PL.HAPCUT2.linked.vcf.gz `ls tmp/hap_*.$tag.PL.HAPCUT2.vcf.gz | sort -k 1,1V`; 
bcftools merge -Oz -o stitch.Herato1603.450000.550000.PL.HAPCUT2.unlinked.vcf.gz `ls tmp/hap_*.$tag.PL.HAPCUT2_unlinked.vcf.gz | sort -k 1,1V`; 
tabix stitch.Herato1603.450000.550000.PL.HAPCUT2.linked.vcf.gz;
tabix stitch.Herato1603.450000.550000.PL.HAPCUT2.unlinked.vcf.gz
perl ../phaseSwitch.pl $tag linked
perl ../phaseSwitch.pl $tag unlinked
cd ..
