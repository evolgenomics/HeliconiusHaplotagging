# Disclaimer
Note that the scripts here is provided here as an actual illustrative example. 
There are many sample name pre- and suffixes that are specific to the dataset. The script as such would not be runnable without modification.

# Note
The key command call for our genotyping is [STITCH](https://github.com/rwdavies/STITCH/). 

```R
STITCH.R --chr=$chr --regionStart=$start --regionEnd=$end --buffer=25000 --bamlist=<BAM_LIST> --sampleNames_file=<SAMPLE_NAME_LIST> --posfile=<POS>\
--outputdir=${out_dir}/${out_stem}_${chr}_${start}_$end --K=30 --method=diploid --nGen=500 --nCores=1 --readAware=TRUE --keepInterimFiles=FALSE \
--shuffle_bin_radius=500 --expRate=5 --iSizeUpperLimit=500000 --keepSampleReadsInRAM=TRUE
```

This particular set of parameters were optimized for Heliconius erato butterflies using known colour genotypes that are segregating as Mendelian loci in this system.
