# Note

This is the code for detecting structural variants

# Basic pipeline

- generate intervals based on the genome assembly, using `bedtools makewindows`, here 10 kbp. See section Interval Writing below.
- populate a pairwise matrix of intervals using the utility code: `interval_write.pl`
- generate a BED file for the samples with linked reads using `bed_write.full.pl`
- index the bgzipped BED file using `tabix`
- on a compute cluster, calculate a Jaccard index for beadTag sharing across windows using the code `struct_jaccard.pl`, looping through the intervals list using an SGE array job, with each job index used to call a joblist file
- basic post-processing of the Jaccord matrix by calculating genome-wide background threshold levels
- analyze the Jaccard matrix to detect any anomolies and calling regions of interest (ROIs)

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
