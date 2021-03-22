#!/bin/bash
tag="example"

echo "***NOTE: The entire pipeline with 968 samples takes a long time to run. Consider running this in a screen or similar virtual terminals.***";
echo "[run_all.sh] Setting up...
";

chmod 755 setup.sh
./setup.sh

echo "[run_mspms.sh] Simulate using mspsm 968 haplotypes, 13505 SNPs, at 10 rho using the tag \"example\"...
";

scripts/run_mspms.sh example 968 13505 10 1

echo "[run_all.sh] Generating truth haploytpes...
";

perl scripts/simulate_haplotypes_ARG.pl $tag 1 968 0.01

echo "[run_all.sh] Starting simulation...
";

scripts/simulate_reads-molecules.sh $tag

scripts/run_STITCH.sh $tag linkedReads simBam/bamfiles.linkedReads.$tag.968.list --use_bx_tag=TRUE  
scripts/run_STITCH.sh $tag shortReads simBam/bamfiles.linkedReads.$tag.968.list

for datatype in linked short; do
    cd STITCH_out/STITCH_968_${datatype}Reads_$tag;
    ../../scripts/post_STITCH.sh $tag STITCH_out/STITCH_968_${datatype}Reads_$tag; 
    ../../scripts/run_hapcut2.sh STITCH_out/STITCH_968_${datatype}Reads_$tag;                           
    cd ..
done;
