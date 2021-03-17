# Simulations
The purpose of this repo is to allow users to simulate our "population haplotyping" pipeline under a range of diversity conditions.
The focal example we use here is from the *Heliconius erato* system. 

Broadly speaking, this simulation scheme creates a set of diploid individuals with distinct haplotypes and allows recombination between them. It then simulates short reads against each haplotype and map these against a canonical reference genome assembly. A separate Perl script would simulate molecules (each marked by the SAM tag "BX" for beadTag) and choose from among the set of mapped reads (as a BAM file) to form linked reads. Then pairs of BX-tagged haploid BAM files would be merged into diploid BAMs. 

At this point the set of BAM files can be optionally subsampled by molecules (BX tag), to reflect varying sequencing coverage.

The entire set of BX-tagged, diploid BAM files (optionally subsampled) would then be placed against the 

For the purpose of this simulation, we will be creating a 1 Mbp region, chosen from a putatively neutral, background location.



# Dependencies
To run this simulation, you'll need:
- [msprime](http://https://github.com/tskit-dev/msprime)
- [dwgsim](https://github.com/nh13/DWGSIM)
- [bwa](https://github.com/lh3/bwa)
- [samtools](https://github.com/samtools/samtools)
- [bcftools](https://github.com/samtools/bcftools)
- [STITCH](https://github.com/rwdavies/STITCH)
- [datamash](https://www.gnu.org/software/datamash/)

Here is a quick overview of the pipelines used in the paper:
![Pipeline](Pipelines_simulations_1.png?raw=true "Simulation pipeline")

