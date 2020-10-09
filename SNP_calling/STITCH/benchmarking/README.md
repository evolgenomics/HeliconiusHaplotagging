# Note

This is the code for testing against Longshanks mice

- The Perl script `compare_concordance.pl` is used to check the per-individual genotype concordance between the STITCH output between the deeper short-read sequences published in Castro et al., eLIFE 2019, and the new haplotag sequences generated in this work.

- In brief, the genotype/dosage GTDS tables from each of 10 subsampled replicates were evaluated against the full sequence, by means of a simple bash loop, like so:

```bash
for gtds in *.GTDS; do\
    echo $gtds"        "`perl compare_concordance.pl $gtds;` \
done
```
# Example output

- Then the output is redirected to a summary file - shown here as an example, with the filename in the first column indicating the run and subsampling conditions, and then in the 4th column, the GT concordance.
# Note on subsampling schemes

- "Original" here refers to using short-read only, without considering linked reads BX-tag information
- "merged" refers to simple subsampling, in which BAM files were subsampled without regard to the BX tag, then merged according to BX tag. This corresponds to the scenario where samples are pooled, but then sequenced at far reduced throughput. This is *not* how multiple haplotagged should be pooled because it will result in far too many molecules/beadTags but with too few reads for each, resulting in many "molecules" with only a single paired-end read.
- "molSub" refers to molecular subsampling, in which BAM files were first sorted by BX tag, and then a subset of BX tags were chosen, retaining all reads carrying the same beadTag. This corresponds to our pooling scheme, in which a subset of beads, hence beadTags, were taken before pooling across samples. This has the effect of retaining the original molecules and achieving a similar molecule/sequencing effort ratio, despite multiplexing and pooling.
