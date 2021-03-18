#!/bin/bash
tag=$1
datatype=$2
bamlist=$3
addPara=$4

if [ ! -e STITCH_out ]; then mkdir STITCH_out; fi
/fml/chones/local/STITCH/STITCH.R --chr=Herato1603 --regionStart=3450000 --regionEnd=3550000 --buffer=25000 --bamlist=$bamlist --sampleNames_file=sourceData/sample.names \
	--posfile=sourceData/PC062_merged_Herato1603.3.45Mb.PL.AD.HAPCUT2.pos \
	--outputdir=STITCH_out/STITCH_968_${datatype}_$tag --K=40 --method=diploid --nGen=500 \
	--nCores=16 --readAware=TRUE --keepInterimFiles=TRUE --shuffle_bin_radius=500 --expRate=5 --iSizeUpperLimit=500000 \
	--keepSampleReadsInRAM=FALSE $addPara 2>&1 |tee log/STITCH_$datatype.$tag.out
