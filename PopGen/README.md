# Introduction
Scripts and files in this folder are example scripts for selection scans and summary statistics between populations.

Briefly, we used: 
- omegaPlus
- sweeD
- selScan
to perform some of the haplotype-based selection scans.

In addition, we used angsd to generate statistics based on site frequency spectra (SFS).

Notes on particular tests are discussed below.

# omegaPlus
Here, we found that setting a small enough grid size is key to recovering the strongest possible signal. This is because the scan for omega_{max} is not exhaustive, therefore at certain coarse grid spacings, the programme can miss very strong omega signals.
Therefore in our pipeline, we dynamically adjusted the grid size based on the total length of the sequence in question.

To represent the phased haplotypes, we opted to use the MaCS format. We therefore added a dummy header in the MaCS.in.template file, and then dynamically recoded the VCF file through `awk` into the MaCS file format.
