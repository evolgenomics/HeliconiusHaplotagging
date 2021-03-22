#!/bin/bash
tag=$1;
folder=$2;

#cd $folder;
for i in `seq 0 2 967`; do echo "HAPCUT - $i"; ../../scripts/hapcutVcf.helera1_demo.sh hap_$i.linkedReads.diploid.$tag.bam stitch.Herato1603.450000.550000.PL.reheadered.vcf.gz $tag; done
bcftools merge -Oz -o stitch.Herato1603.450000.550000.PL.HAPCUT2.linked.vcf.gz `ls tmp/hap_*.$tag.PL.HAPCUT2.vcf.gz | sort -k 1,1V`; 
bcftools merge -Oz -o stitch.Herato1603.450000.550000.PL.HAPCUT2.unlinked.vcf.gz `ls tmp/hap_*.$tag.PL.HAPCUT2_unlinked.vcf.gz | sort -k 1,1V`; 
tabix stitch.Herato1603.450000.550000.PL.HAPCUT2.linked.vcf.gz;
tabix stitch.Herato1603.450000.550000.PL.HAPCUT2.unlinked.vcf.gz
perl ../../scripts/phaseSwitch_highConf_quant.pl $tag linked
perl ../../scripts/phaseSwitch_highConf_quant.pl $tag unlinked
#cd ..
