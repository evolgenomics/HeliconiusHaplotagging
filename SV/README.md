# Note

This is the code for detecting structural variants

# Basic pipeline

- generate intervals based on the genome assembly, using `bedtools makewindows`, here 10 kbp. See section Interval Writing below.
- populate a pairwise matrix of intervals using the utility code: `interval_write.pl`
- generate a BED file for the samples with linked reads using `bed_write.full.pl`
- index the bgzipped BED file using `tabix`
- on a compute cluster, calculate a Jaccard index for beadTag sharing across windows using the code `struct_jaccard.pl`, looping through the intervals list using an SGE array job, with each job index used to call a joblist file
- basic post-processing of the Jaccord matrix by calculating genome-wide background threshold levels
- analyze the Jaccard matrix to detect any anomolies and calling regions of interest (ROIs), here using `pick_ROI.pl`

# Interval writing

A basic example is shown here on how to write out the per-chromosome pairwise interval job list:
```bash
#Writing out the intervals for Hmel2.5 - making a per chromosome interval file, in the CHROM:START-END format. Here the first 7 characters of col1 is the hmel2.5 scaffold name.
awk '{print $1":"$2"-"$3 > "hmel2_5.scaffold.10kb_win."substr($1,1,7)".intervals"}'  hmel2_5.scaffold.10kb_win.bed

#Concatenating also the unassembled scaffolds onto each file with a bash loop
for i in `seq 21`; do pi=`printf %02d $i`; cat hmel2_5.scaffold.10kb_win.Hmel200.intervals >> hmel2_5.scaffold.10kb_win.Hmel2$pi.intervals; done

#Then running the utility perl script interval_write.pl to create the bulk joblist file
for i in hmel2_5.scaffold.10kb_win.Hmel*.intervals; do perl interval_write.pl $i; done > joblist/hmel2.5.bulk

#Then we split the bulk joblist files into chunks, by line.
```

# Calling ROIs
The ROIs are determined by the extent of beadTag sharing across two windows, along with some information from a genome-wide normalized background sharing rate. 

The background rate is calculated this way, using helera1_demo as an example:
```bash
#First adding a column to the jaccard_matrix output file by expression the fraction of shared BX tags over all BX tag combinations - c12 = c11/(c10*c9)
zcat hmel2.5.10k_win.Hmel2*.sites.jaccard_matrix.gz | cut -f 1-11 | awk 'BEGIN {OFS="\t"}; {$12 = $11/($10*$9); print $0}' | datamash median 12
##2.57147e-08

#This value is stored into the $background_rate variable
background_rate=`zcat hmel2.5.10k_win.Hmel2*.sites.jaccard_matrix.gz  | cut -f 1-11 | awk 'BEGIN {OFS="\t"}; {$12 = $11/($10*$9); print $0}'  | datamash median 12`

#Now make the table
zcat hmel2.5.10k_win.Hmel2*.sites.jaccard_matrix.gz | cut -f 1-11 | awk 'BEGIN {OFS="\t"}; {$12 = $11/($10*$9);$13=$12/'$background_rate'; print $0}' > hmel2.5.10k_win.Hmel2*.sites.jaccard_matrix.front
```

# Misc
Accumulative offset file
```bash
awk '{print substr($1,1,8)"\t"$0}' helera1_demo.accum.sizes | datamash groupby 1 first 5 > helera1_demo.accum.offset
```
