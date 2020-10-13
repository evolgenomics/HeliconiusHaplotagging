# Note

This is the code for detecting structural variants

# Basic pipeline

- generate intervals based on the genome assembly, using `bedtools makewindows`
- populate a pairwise matrix of intervals using the utility code: `interval_write.pl`
- generate a BED file for the samples with linked reads using `bed_write.full.pl`
- index the bgzipped BED file using `tabix`
- on a compute cluster, calculate a Jaccard index for beadTag sharing across windows using the code `struct_jaccard.pl`, looping through the intervals list using an SGE array job, with each job index used to call a joblist file
- basic post-processing of the Jaccord matrix by calculating genome-wide background threshold levels
- analyze the Jaccard matrix to detect any anomolies and calling regions of interest (ROIs)
