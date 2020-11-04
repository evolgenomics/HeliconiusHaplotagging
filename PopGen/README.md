# Introduction
Scripts and files in this folder are example scripts for selection scans and summary statistics between populations.

Briefly, we used: 
- omegaPlus
- sweeD
- selScan
to perform some of the haplotype-based selection scans.

We computed pi with vcftools for delta pi.

In addition, we used angsd to generate statistics based on site frequency spectra (SFS) such as FST. 

Lastly, we computed population structure with NGSadmix and performed genomewide association mapping (GWAS) with angsd.

Notes on particular tests are discussed below.

# FST
We computed FST per-site and in 10 kb windows with angsd using the genotype likelihoods emitted by HAPCUT2. We then used a Hidden Markov Model (HMM) approach to infer regions of particularly high FST indicative of divergent selection.

# Delta pi
Here, we computed the nucleotide diversity (pi) for each race separately with vcftools. Whereas pi can vary due to differences in mutation rates, recombination rates, etc, these factors should affect both races similarly. However, selection in only one race is expected to lead to a decreased pi in that race and thus a difference in pi. Delta pi, i.e. the difference in pi between the two races compared is thus expected to deviate from zero in genomic regions where one race experienced more selection than the other one.

# omegaPlus
Here, we found that setting a small enough grid size is key to recovering the strongest possible signal. This is because the scan for omega_{max} is not exhaustive, therefore at certain coarse grid spacings, the programme can miss very strong omega signals.
Therefore in our pipeline, we dynamically adjusted the grid size based on the total length of the sequence in question.

To represent the phased haplotypes, we opted to use the MaCS format. We therefore added a dummy header in the MaCS.in.template file, and then dynamically recoded the VCF file through `awk` into the MaCS file format.

