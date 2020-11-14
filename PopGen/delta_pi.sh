#!/bin/bash

module load htslib/1.2.1 bcftools-1.9-gcc-5.4.0-b2hdt5n samtools/1.3.1 vcftools

# Compute 10 kb pi for each scaffold
for file in ../vcf/run164_merged_*.PL.AD.HAPCUT2.inf0.05.vcf.gz
do 
  prefix=`basename $i`
  prefix=${prefix%.vcf.gz}
  vcftools --gzvcf $file --keep malleti.samples --window-pi 10000 --window-pi-step 2500 --out $prefix.malleti.10kb
  vcftools --gzvcf $file --keep plesseni.samples --window-pi 10000 --window-pi-step 2500 --out $prefix.plesseni.10kb
done

# Compute 50 kb pi for each scaffold
for file in ../vcf/run164_merged_*.PL.AD.HAPCUT2.inf0.05.vcf.gz
do 
  prefix=`basename $i`
  prefix=${prefix%.vcf.gz}
  vcftools --gzvcf $file --keep malleti.samples --window-pi 50000 --out $prefix.malleti.50kb
  vcftools --gzvcf $file --keep plesseni.samples --window-pi 50000 --out $prefix.plesseni.50kb
done

# Delta pi was then computed in R as pi(highland subspecies) - pi(lowland subspecies)
# For parts of Fig3, see in https://github.com/evolgenomics/HeliconiusHaplotagging/blob/main/PopGen/plotGenomeScans_Fig3.r

# R code for Supplementary figure showing how delta pi of the four colour loci compares to the genome-wide distribution

R

# Read in the datasets for all four species:
thetaLativ50kb<-read.table("PC062_merged.PL.AD.HAPCUT2.inf0.05.lativitta.50kb.windowed.pi",header=T)
thetaNota50kb<-read.table("PC062_merged.PL.AD.HAPCUT2.inf0.05.notabilis.50kb.windowed.pi",header=T)
thetaMalle50kb<-read.table("run164_merged.PL.AD.HAPCUT2.inf0.05.malleti.50kb.windowed.pi",header=T)
thetaPless50kb<-read.table("run164_merged.PL.AD.HAPCUT2.inf0.05.plesseni.50kb.windowed.pi",header=T)

# Make the header more readable
names(thetaLativ50kb)<-c("Chr","start","end","n","pi")
names(thetaNota50kb)<-c("Chr","start","end","n","pi")
names(thetaMalle50kb)<-c("Chr","start","end","n","pi")
names(thetaPless50kb)<-c("Chr","start","end","n","pi")

# add window centre coordinates
thetaLativ50kb$WinCenter<-thetaLativ50kb$start+25000
thetaNota50kb$WinCenter<-thetaNota50kb$start+25000
thetaMalle50kb$WinCenter<-thetaMalle50kb$start+25000
thetaPless50kb$WinCenter<-thetaPless50kb$start+25000

# Extract WntA scaffold
thetaLativ50kbWntA<-thetaLativ50kb[thetaLativ50kb$Chr=="Herato1001",]
thetaNota50kbWntA<-thetaNota50kb[thetaNota50kb$Chr=="Herato1001",]
thetaMalle50kbWntA<-thetaMalle50kb[thetaMalle50kb$Chr=="Hmel210001o",]
thetaPless50kbWntA<-thetaPless50kb[thetaPless50kb$Chr=="Hmel210001o",]

# Extract vvl scaffold
thetaLativRo50kb<-thetaLativ50kb[thetaLativ50kb$Chr=="Herato1301",]
thetaNotaRo50kb<-thetaNota50kb[thetaNota50kb$Chr=="Herato1301",]
thetaMalle50kbRo<-thetaMalle50kb[thetaMalle50kb$Chr=="Hmel213001o",]
thetaPless50kbRo<-thetaPless50kb[thetaPless50kb$Chr=="Hmel213001o",]

# Extract optix scaffold
thetaLativ50kbOptix<-thetaLativ50kb[thetaLativ50kb$Chr=="Herato1801",]
thetaNota50kbOptix<-thetaNota50kb[thetaNota50kb$Chr=="Herato1801",]
thetaMalle50kbOptix<-thetaMalle50kb[thetaMalle50kb$Chr=="Hmel218003o",]
thetaPless50kbOptix<-thetaPless50kb[thetaPless50kb$Chr=="Hmel218003o",]

# Extract cortex scaffold
thetaLativ50kbN<-thetaLativ50kb[thetaLativ50kb$Chr=="Herato1505",]
thetaNota50kbN<-thetaNota50kb[thetaNota50kb$Chr=="Herato1505",]
thetaMalle50kbN<-thetaMalle50kb[thetaMalle50kb$Chr=="Hmel215003o",]
thetaPless50kbN<-thetaPless50kb[thetaPless50kb$Chr=="Hmel215003o",]



# Get the window with the lowest delta pi in the regions of optix, WntA, vvl and cortex
# The region ranges are those shown in Fig. 3b and 3c
# They are inferred here: https://github.com/evolgenomics/HeliconiusHaplotagging/blob/main/PopGen/plotGenomeScans_Fig3.r

# optix melpomene
rangeMelOptix<-c(0,1787500)
optixMel<-min(thetaPless50kb[thetaPless50kbOptix$start>rangeMelOptix[1]&thetaPless50kbOptix$end<rangeMelOptix[2],"pi"]-
                thetaMalle50kb[thetaMalle50kbOptix$start>rangeMelOptix[1]&thetaMalle50kbOptix$end<rangeMelOptix[2],"pi"])
# WntA melpomene
rangeWntA<-c(2335500,4335500)
WntAMel<-min(thetaPless50kb[thetaPless50kbWntA$start>rangeMelWntA[1]&thetaPless50kbWntA$end<rangeMelWntA[2],"pi"]-
               thetaMalle50kb[thetaMalle50kbWntA$start>rangeMelWntA[1]&thetaMalle50kbWntA$end<rangeMelWntA[2],"pi"])
# vvl melpomene
rangeMelRo<-c(9352500,11352500)
RoMel<-min(thetaPless50kb[thetaPless50kbRo$start>rangeMelRo[1]&thetaPless50kbRo$end<rangeMelRo[2],"pi"]-
             thetaMalle50kb[thetaMalle50kbRo$start>rangeMelRo[1]&thetaMalle50kbRo$end<rangeMelRo[2],"pi"])
# cortex melpomene
rangeMelN<-c(442500,2442500)
NMel<-min(thetaPless50kbN[thetaPless50kbN$start>rangeMelN[1]&thetaPless50kbN$end<rangeMelN[2],"pi"]-
            thetaMalle50kbN[thetaMalle50kbN$start>rangeMelN[1]&thetaMalle50kbN$end<rangeMelN[2],"pi"])
# optix erato
rangeEraOptix<-c(385000,2385000)
optixEra<-min(thetaNota50kbOptix[thetaNota50kbOptix$start>rangeEraOptix[1]&thetaNota50kbOptix$end<rangeEraOptix[2],"pi"]-
              thetaLativ50kbOptix[thetaLativ50kbOptix$start>rangeEraOptix[1]&thetaLativ50kbOptix$end<rangeEraOptix[2],"pi"])
# WntA erato
rangeEraWntA<-c(3645000,5645000)
WntAEra<-min(thetaNota50kbWntA[thetaNota50kbWntA$start>rangeEraWntA[1]&thetaNota50kbWntA$end<rangeEraWntA[2],"pi"]-
              thetaLativ50kbWntA[thetaLativ50kbWntA$start>rangeEraWntA[1]&thetaLativ50kbWntA$end<rangeEraWntA[2],"pi"])
# vvl erato
rangeEraRo<-c(13345000,15345000)
RoEra<-min(thetaNotaRo50kb[thetaNotaRo50kb$start>rangeEraRo[1]&thetaNotaRo50kb$end<rangeEraRo[2],"pi"]-
              thetaLativRo50kb[thetaLativRo50kb$start>rangeEraRo[1]&thetaLativRo50kb$end<rangeEraRo[2],"pi"])
# cortex erato
rangeEraN<-c(1115000,3115000)
NEra<-min(thetaNota50kbN[thetaNota50kbN$start>rangeEraN[1]&thetaNota50kbN$end<rangeEraN[2],"pi"]-
              thetaLativ50kbN[thetaLativ50kbN$start>rangeEraN[1]&thetaLativ50kbN$end<rangeEraN[2],"pi"])


# Compute one-sided empirical p-value as the proportion of more extreme values than the one at the colour locus
require(scales)

# H. melpomene
# optix
pOpM<-scientific(length(which((thetaPless50kb$pi-thetaMalle50kb$pi)<optixMel))/length(thetaMalle50kb$pi))
# WntA
pWntM<-scientific(length(which((thetaPless50kb$pi-thetaMalle50kb$pi)<WntAMel))/length(thetaMalle50kb$pi))
# vvl
pRoM<-scientific(length(which((thetaPless50kb$pi-thetaMalle50kb$pi)<RoMel))/length(thetaMalle50kb$pi))
# cortex
pNM<-scientific(length(which((thetaPless50kb$pi-thetaMalle50kb$pi)<NMel))/length(thetaMalle50kb$pi))

# H. erato
# optix
pOpE<-scientific(length(which((thetaNota50kb$pi-thetaLativ50kb$pi)<optixEra))/length(thetaLativ50kb$pi))
# WntA
pWntE<-scientific(length(which((thetaNota50kb$pi-thetaLativ50kb$pi)<WntAEra))/length(thetaLativ50kb$pi))
# vvl
pRoE<-scientific(length(which((thetaNota50kb$pi-thetaLativ50kb$pi)<RoEra))/length(thetaLativ50kb$pi))
# cortex
pNE<-scientific(length(which((thetaNota50kb$pi-thetaLativ50kb$pi)<NEra))/length(thetaLativ50kb$pi))


# Plot delta pi distribution

par(mfrow=c(2,1),mgp=c(2.2,0.5,0),mar=c(4,3.2,1,1))
# Genome-wide distribution
plot(density(thetaNota50kb$pi-thetaLativ50kb$pi),xlim=c(-0.017,0.017),lwd=2,
     xlab=expression(Delta~pi~"("~pi~"notabilis - "~pi~"lativitta)"),main="")
# Add vertical lines for the colour loci
abline(v=c(optixEra,WntAEra,RoEra,NEra),col=c("#C0202A","#CCCC00","#0000CC","#CC33FF"),lwd=2)
# Add a legend also showing the p-values
legend("topright",legend=c(paste0("optix p=",pOpE),paste0("WntA p=",pWntE),
                           paste0("Ro p=",pRoE),paste0("cortex p=",pNE)),
       col=c("#C0202A","#CCCC00","#0000CC","#CC33FF"),lwd=2,bty="n")

# Genome-wide distribution
plot(density(thetaPless50kb$pi-thetaMalle50kb$pi),xlim=c(-0.01,0.01),lwd=2,
     xlab=expression(Delta~pi~"("~pi~"plesseni - "~pi~"malleti)"),main="")
# Add vertical lines for the colour loci
abline(v=c(optixMel,WntAMel,RoMel,NMel),col=c("#C0202A","#CCCC00","#0000CC","#CC33FF"),lwd=2)
# Add a legend also showing the p-values
legend("topright",legend=c(paste0("optix p=",pOpM),paste0("WntA p=",pWntM),
                           paste0("Ro p=",pRoM),paste0("cortex p=",pNM)),
       col=c("#C0202A","#CCCC00","#0000CC","#CC33FF"),lwd=2,bty="n")




