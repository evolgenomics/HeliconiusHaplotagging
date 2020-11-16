# Introduction
Scripts and files in this folder are example scripts for selection scans, trait association mapping and summary statistics between populations.

Briefly, we used: 
* angsd for [FST](https://github.com/evolgenomics/HeliconiusHaplotagging/blob/main/PopGen/FST_melpomene.sh)
* vcftools for [delta pi](https://github.com/evolgenomics/HeliconiusHaplotagging/blob/main/PopGen/delta_pi.sh)

and for haplotype-based selection scans:
* omegaPlus
* sweeD
* selScan

Genome-wide association mapping (GWAS) and interval mapping
* [population structure](https://github.com/evolgenomics/HeliconiusHaplotagging/blob/main/PopGen/NGSadmix.sh) used as a covariate for GWAS was inferred with NGSadmix 
* [GWAS](https://github.com/evolgenomics/HeliconiusHaplotagging/blob/main/PopGen/GWAS) was performed with angsd 
* interval mapping was performed using a custom R [script](https://github.com/evolgenomics/HeliconiusHaplotagging/blob/main/PopGen/interval_mapping_melpomene.sh)

The R script for Fig. 3 combining these statistics is given [here](https://github.com/evolgenomics/HeliconiusHaplotagging/blob/main/PopGen/plotGenomeScans_Fig3.r).

Notes on particular tests are given below.

# FST
We computed FST per-site and in 10 kb windows with angsd using the genotype likelihoods emitted by HAPCUT2. We then used a Hidden Markov Model (HMM) approach to infer regions of particularly high FST indicative of divergent selection.

# Delta pi
Here, we computed the nucleotide diversity (pi) for each race separately with vcftools. Whereas pi can vary due to differences in mutation rates, recombination rates, etc, these factors should affect both races similarly. However, selection in only one race is expected to lead to a decreased pi in that race and thus a difference in pi. Delta pi, i.e. the difference in pi between the two races compared is thus expected to deviate from zero in genomic regions where one race experienced more selection than the other one.

# omegaPlus
Here, we found that setting a small enough grid size is key to recovering the strongest possible signal. This is because the scan for omega_{max} is not exhaustive, therefore at certain coarse grid spacings, the programme can miss very strong omega signals.
Therefore in our pipeline, we dynamically adjusted the grid size based on the total length of the sequence in question.

To represent the phased haplotypes, we opted to use the MaCS format. We therefore added a dummy header in the MaCS.in.template file, and then dynamically recoded the VCF file through `awk` into the MaCS file format.

# GWAS
We inferred the population structure with [NGSadmix](http://www.popgen.dk/software/index.php/NgsAdmix) and used the admixture proportions as covariates for the GWAS (genome-wide association study) in [angsd](http://www.popgen.dk/angsd/index.php/Association). The optimal number of K clusters was inferred with [Clumpak](http://clumpak.tau.ac.il/). We performed GWAS using the Score test in angsd for each phenotypic trait separately.

# Interval mapping
To refine specific regulatory elements or tightly linked genes in Heliconius melpomene, we performed interval mapping for traits that mapped to optix and cortex and/or nearby genes. This allowed us to narrow down which regulatory element of optix controls which aspect of the red patterns and it allowed us to figure out that the variable presence of a yellow spot at the forewing base maps to cortex, whereas the distribution of red scales maps to two overlapping genes (domeless/washout) that are tighly linked to cortex.
