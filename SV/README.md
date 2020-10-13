# Note

This is the code for detecting structural variants

# Basic pipeline

- generate intervals based on the genome assembly, using `bedtools makewindows`, here 10 kbp. See section Interval Writing below.
- populate a pairwise matrix of intervals using the utility code: `interval_write.pl`
- generate a BED file for the samples with linked reads using `bed_write.full.pl`
- index the bgzipped BED file using `tabix`

Continguity SV (Inversions/Transpositions/Translocations):
- on a compute cluster, calculate a Jaccard index for beadTag sharing across windows using the code `struct_jaccard.pl`, looping through the intervals list using an SGE array job, with each job index used to call a joblist file
- basic post-processing of the Jaccord matrix by calculating genome-wide background threshold levels
- analyze the Jaccard matrix to detect any anomolies and calling regions of interest (ROIs), here using `pick_ROI.pl`

Size SV (Insertions/deletions):
For indels we only need to handle a single interval. Strictly speaking indels should appear as anomalies in both size *and* beadTag sharing. This section complements the information by analyzing the molecule size
- here we use a different cluster script `struct_indel.sge.sh`. The main difference is that it calls `structural_indel.pl` instead.
- These results are then aggregated by a simple bash loop:
```bash
for i in `seq 727`; do pi=`printf %04d $i`; cat helera1_demo.10k_win.interval_$pi.molSize.out; done |  awk '{print $0 > "helera1_demo.10k_win."substr($2,1,8)".molSize.out"}'
#Add header
for i in `seq 21`; do pi=`printf %02d $i`; head -1 molSize.header  | xargs -i sed -i '1s;^;'{}'\n;' helera1_demo.10k_win.Herato$pi.molSize.out ; done
```
- Then the values are normalized per population group, and used to create a reformatted normalized table:
```bash
#Genome-wide average molecule size per population:
awk '!/INTERVAL/' helera1_demo.10k_win.Herato*.molSize.out | datamash mean `seq 7 3 67 | paste -d"," -s` q1 `seq 7 3 67 | paste -d"," -s` q3 `seq 7 3 67 | paste -d"," -s` --narm;
 
##9675.1512597888 8824.4990473382 6696.8443347322 13585.936545682 10002.271430342 13110.006695393 10633.376510753 9726.3031480578 11057.949014198 12722.594348395 10801.981998112 17603.327013595 11951.680434555 10716.22759819 11899.013809093 9577.0992904137 6968.1988056376 7463.3213813809 9149.130727352 8751.0966154174 8245.6221745849

#Then reform the header
for i in `seq 21`; do pi=`printf %02d $i`; awk 'BEGIN {split("9675.1512597888 8824.4990473382 6696.8443347322 13585.936545682 10002.271430342 13110.006695393 10633.376510753 9726.3031480578 11057.949014198 12722.594348395 10801.981998112 17603.327013595 11951.680434555 10716.22759819 11899.013809093 9577.0992904137 6968.1988056376 7463.3213813809 9149.130727352 8751.0966154174 8245.6221745849",mean,/ +/);OFS="\t";}; !/INTERVAL/ {for(i=7;i<=NF;i+=3){$i=$i"\t"$i/(mean[(i-4)/3])};print $0}' helera1_demo.10k_win.Herato$pi.molSize.out; done | sort -k 2,2 -k 3,3V | uniq > helera1_demo.10k_win.all.molSize.normalized.out
cat molSize.header.normalized | xargs -i sed -i '1s;^;'{}'\n;' helera1_demo.10k_win.all.molSize.normalized.out

#Clean-up 
rm helera1_demo.10k_win.interval_????.molSize.out

#Making a list of importnt internal indels in the genome - here, by picking out the first and last 50 kbp of a scaffold, and putting these into a set of censored grep patterns
 grep ":1-" -A 4 ../helera1_demo.10k_win.Herato*.intervals | cut -c 46- | awk '/Herato/' > borders.grep
 grep -v "0000$" -B 4 ../helera1_demo.10k_win.Herato*.1.intervals | cut -c 46- | awk '/Herato/' >> borders.grep 

#The pluck out the columns that corresponds to all samples (c65), or just notabilis (c66 from Site 1) or lativitta (c82). Then normalize these against the genome-wide averages of these columns, and also producing a "delta" column. These positions are filtered against the borders.grep list.
#The final intervals that pass the 5%/95% quantile thresholds: 0.700 and 2.027 are then dumped out to the Indel_internal_global.list
cut -f 1,65,66,82,83-86 helera1_demo.10k_win.Herato*.indel | sort -k 2,2nr | awk '{print $1"\t"$2/22.950946072831"\t"$3/1.3127928815398"\t"$4/0.52134928880377"\t"$3/1.3127928815398-$4/0.52134928880377}' | sort -k 2,2n | grep -v -f borders.grep | awk '$2 <= 0.700683; $2 >= 2.0271' > Indel_internal_global.list

#This is then sorted by scaffold and positions, and clustered to create a set of cohoerent loci positions [415 for H. erato] that contain putative indels.
sed 's/:/\t/;s/-/\t/g' Indel_internal_global.list | sort -k 1,1 -k 2,2n | cut -f 1-4  | bedtools cluster -i - -d 50000 | datamash groupby 5 first 1 min 2 max 3 collapse 4 | cut -f 2- > Indel_internal_global.clustered.list 
```

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

The background rate is calculated this way:
```bash
#First adding a column to the jaccard_matrix output file by expression the fraction of shared BX tags over all BX tag combinations - c12 = c11/(c10*c9)
zcat hmel2.5.10k_win.Hmel2*.sites.jaccard_matrix.gz | cut -f 1-11 | awk 'BEGIN {OFS="\t"}; {$12 = $11/($10*$9); print $0}' | datamash median 12
##2.57147e-08

#This value is stored into the $background_rate variable
background_rate=`zcat hmel2.5.10k_win.Hmel2*.sites.jaccard_matrix.gz  | cut -f 1-11 | awk 'BEGIN {OFS="\t"}; {$12 = $11/($10*$9); print $0}'  | datamash median 12`
```

Now make the table
```bash
zcat hmel2.5.10k_win.Hmel2*.sites.jaccard_matrix.gz | cut -f 1-11 | awk 'BEGIN {OFS="\t"}; {$12 = $11/($10*$9);$13=$12/'$background_rate'; print $0}' > hmel2.5.10k_win.Hmel2*.sites.jaccard_matrix.front
```

# Misc
Accumulative offset file
```bash
awk '{print substr($1,1,8)"\t"$0}' helera1_demo.accum.sizes | datamash groupby 1 first 5 > helera1_demo.accum.offset
```
