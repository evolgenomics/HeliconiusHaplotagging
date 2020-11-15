#!/bin/bash 

# Code for interval mapping of the optix region and the cortex region

##################################################################################################
# Get optix region genotypes
##################################################################################################


# Filter the genotypes more stringently by adding a Genotype Quality tag:
# Add GQ tags (=second smallest PL value, capped at 99)

zcat $file.vcf.gz | \
awk  -v OFS='\t' '{if($1~"#"){print $0}
 else{$9=$9":GQ";
  for (i = 10; i <= NF; ++i){
    split($i,a,":");split(a[5],pl,",");
    if(pl[1]>pl[2]){gq=pl[1]; min=2}else{gq=pl[2];min=1} 
    if(pl[3]<=gq){if(pl[3]>pl[min]) gq=pl[3]; else gq=pl[min]};
    if(gg>99) gq=99;
    $i=$i":"gq; 
  } print $0
 }
}' | gzip > $file.GQ.vcf.gz

# Filter out genotypes with GQ<15
vcftools --gzvcf $file.GQ.vcf.gz --minGQ 15 --recode \
--stdout --max-missing 0.75 | gzip > $file.minGQ15.max0.25N.vcf.gz

# Get genotypes counts for the samples  of the different optix phenotypic groups
vcftools --gzvcf $vcfFile --keep optixM.samples --out optix.m --hardy
vcftools --gzvcf $vcfFile --keep optixP.samples --out optix.p --hardy
vcftools --gzvcf $vcfFile --keep optixHet.samples --out optix.het --hardy
vcftools --gzvcf $vcfFile --keep weirdo.samples --out optix.weirdo --hardy  # note, these are the individuals with dennis but no rays

# Concatenate the files of genotype counts for each group
paste <(cut -f 1-3 optix.m.hwe) <(cut -f 3 optix.p.hwe) <(cut -f 3 optix.het.hwe ) <(cut -f 3 optix.weirdo.hwe) > optix.genotypes

# Fix the header so that the group names are represented
sed -i -e 's#OBS(HOM1/HET/HOM2)\tOBS(HOM1/HET/HOM2)\tOBS(HOM1/HET/HOM2)\tOBS(HOM1/HET/HOM2)#HOM1_L/HET_L/HOM2_L\tHOM1_N/HET_N/HOM2_N\tHOM1_het/HET_het/HOM2_het\tHOM1_w/HET_w/HOM2_w#g' \
 -e 's#/#\t#g' optix.genotypes

# Add proportion of HOM1 counts among all genotype counts
awk 'BEGIN{print "HOM1latProp\tHOM1notProp\tHEThetProp"}{if(($6+$7+$8)>0) print $3/($3+$4+$5),$6/($6+$7+$8),$10/($9+$10+$11)}' optix.genotypes > propHom1
paste optix.genotypes propHom1 > optix.genotypes.prop

# Find "diagnostic" positions with allele frequency difference of at least 0.9 between the malleti-like and the plesseni-like groups
awk '{if(($15-$16>0.9 || $16-$15>0.9)&&($9<=10 && $11<=10)) print $0}' optix.genotypes.prop > optix.diagnostic.SNPs.0.9hom.max10het

# Extract those positions from the vcf file
suffix=diagnostic.SNPs.0.9hom.max10het
vcftools --keep optixHet.samples --recode --positions optix.$suffix --gzvcf $vcfFile --out optixhet.$suffix
vcftools --keep optixM.samples --recode --positions optix.$suffix --gzvcf $vcfFile --out optixM.$suffix
vcftools --keep optixP.samples --recode --positions optix.$suffix --gzvcf $vcfFile --out optixP.$suffix
vcftools --keep weirdo.samples --recode --positions optix.$suffix --gzvcf $vcfFile --out optixW.$suffix # note, these are the individuals with dennis but no rays

# Extract the genotypes as 0/0, 0/1, 1/1 

grep ^#CH optixM.$suffix.recode.vcf | cut -f 1,2,4,5,10-1000 > optixM.vcf
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' optixM.$suffix.recode.vcf >> optixM.vcf

grep ^#CH optixP.$suffix.recode.vcf | cut -f 1,2,4,5,10-1000 > optixP.vcf
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' optixP.$suffix.recode.vcf >> optixP.vcf

grep ^#CH optixhet.$suffix.recode.vcf | cut -f 1,2,4,5,10-1000 > optixHet.vcf
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' optixhet.$suffix.recode.vcf >> optixHet.vcf

grep ^#CH optixW.$suffix.recode.vcf | cut -f 1,2,4,5,10-1000 > optixW.vcf
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' optixW.$suffix.recode.vcf >> optixW.vcf

# Combine all vcf files
paste optixM.vcf <(cut -f 5-1000 optixHet.vcf) <(cut -f 5-1000 optixP.vcf) <(cut -f 5-1000 optixW.vcf) > optix.$suffix.vcf

# Transform genotypes to single-digit codes (0|0=0, 0|1=1, 1|1=2)
sed -i -e 's#/#|#g' optix.$suffix.vcf
sed -i -e 's/0|0/0/g' -e 's/0|1/1/g' -e 's/1|0/1/g' -e 's/1|1/2/g' -e 's/.|./NA/g' -e 's/#//' optix.$suffix.vcf


##################################################################################################
# Get cortex region genotypes
##################################################################################################

# Get a vcf file of scaffold 15003o with high quality sites and genotypes only:
bcftools view -i "INFO_SCORE >= 0.7"  -Oz -o sc15003o.inf0.7.vcf.gz \
../vcf/run164_merged_Hmel215003o.PL.AD.HAPCUT2.inf0.05.corr.vcf.gz

# Add GQ tags (=second smallest PL value, capped at 99)
zcat sc15003o.inf0.7.vcf.gz | \
awk  -v OFS='\t' '{if($1~"#"){print $0}
 else{$9=$9":GQ";
  for (i = 10; i <= NF; ++i){
    split($i,a,":");split(a[5],pl,",");
    if(pl[1]>pl[2]){gq=pl[1]; min=2}else{gq=pl[2];min=1} 
    if(pl[3]<=gq){if(pl[3]>pl[min]) gq=pl[3]; else gq=pl[min]};
    if(gg>99) gq=99;
    $i=$i":"gq; 
  } print $0
 }
}' | gzip > sc15003o.inf0.7.GQ.vcf.gz

# Filter out bad genotypes and sites with high missing data proportion
vcftools --gzvcf sc15003o.inf0.7.GQ.vcf.gz --minGQ 15 --max-missing 0.75 \
--recode --stdout | gzip > sc15003o.inf0.7.minGQ15.max0.25N.vcf.gz

vcfFile="sc15003o.inf0.7.minGQ15.max0.25N.vcf.gz" # more stringently filtered. 

# Get genotypes counts of individuals that are fully highland/plesseni-like (HH) or lowland/malleti-like (LL)
vcftools --gzvcf $vcfFile --keep cortex.HH --out cortex.H --hardy
vcftools --gzvcf $vcfFile --keep cortex.LL --out cortex.L --hardy

# Combine the two files of genotype counts
paste <(cut -f 1-3 cortex.L.hwe) <(cut -f 3 cortex.H.hwe) > cortex.genotypes
sed -i -e 's#OBS(HOM1/HET/HOM2)\tOBS(HOM1/HET/HOM2)#HOM1_L/HET_L/HOM2_L\tHOM1_H/HET_H/HOM2_H#g' \
 -e 's#/#\t#g' cortex.genotypes

# Add the proportion of HOM1 genotypes
awk 'BEGIN{print "HOM1_H_Prop\tHOM1_L_Prop"}{if(($6+$7+$8)>0) print $3/($3+$4+$5),$6/($6+$7+$8)}' cortex.genotypes > propHom1
paste cortex.genotypes propHom1 > cortex.genotypes.prop

# Find "diagnostic" positions (allele frequency difference is at least 0.9)
awk '{if(($9-$10>=0.9 || $10-$9>=0.9)&&($4<=5 && $7<=5)) print $0}' cortex.genotypes.prop > cortex.diagnostic.SNPs.0.9hom.max5het#suffix=diagnostic.SNPs

suffix=diagnostic.SNPs.0.9hom.max5het

# Extract those positions from the vcf file for each phenotypic group 
# (e.g. cortex.HH are individuals with highland like absence of spot and highland-like red distribution)
vcftools --keep cortex.HL --recode --positions cortex.$suffix --gzvcf $vcfFile --out cortexHL.$suffix
vcftools --keep cortex.LH --recode --positions cortex.$suffix --gzvcf $vcfFile --out cortexLH.$suffix
vcftools --keep cortex.HH --recode --positions cortex.$suffix --gzvcf $vcfFile --out cortexHH.$suffix
vcftools --keep cortex.LL --recode --positions cortex.$suffix --gzvcf $vcfFile --out cortexLL.$suffix

# Strip away all non-necessary information and show only the genotypes as 0/0, 0/1, 1/1

grep ^#CH cortexLL.$suffix.recode.vcf | cut -f 1,2,4,5,10-1000 > cortexLL.vcf
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' cortexLL.$suffix.recode.vcf >> cortexLL.vcf

grep ^#CH cortexHH.$suffix.recode.vcf | cut -f 1,2,4,5,10-1000 > cortexHH.vcf
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' cortexHH.$suffix.recode.vcf >> cortexHH.vcf

grep ^#CH cortexLH.$suffix.recode.vcf | cut -f 1,2,4,5,10-1000 > cortexLH.vcf
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' cortexLH.$suffix.recode.vcf >> cortexLH.vcf

grep ^#CH cortexHL.$suffix.recode.vcf | cut -f 1,2,4,5,10-1000 > cortexHL.vcf
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' cortexHL.$suffix.recode.vcf >> cortexHL.vcf

paste cortexHH.vcf <(cut -f 5-1000 cortexHL.vcf) <(cut -f 5-1000 cortexLL.vcf) <(cut -f 5-1000 cortexLH.vcf) > cortex.$suffix.vcf

# replace the genotypes by single-digit value 
sed -i -e 's#/#|#g' cortex.$suffix.vcf
sed -i -e 's/0|0/0/g' -e 's/0|1/1/g' -e 's/1|0/1/g' -e 's/1|1/2/g' -e 's/.|./NA/g' -e 's/#//' cortex.$suffix.vcf


##################################################################################################

# R code to plot haplotypes for interval mapping of the optix and the cortex regions

##################################################################################################

# optix

# Read in the genotypes at all sites with an allele frequency difference of 0.9 between individuals of highland and lowland optix phenotypes
genotypes<-read.table("optix.diagnostic.SNPs.0.9hom.max10het.vcf",header=T,sep="\t")


# Specify plesseni allele as 1
genotypes[rowMeans(genotypes[,5:104],na.rm=T)>rowMeans(genotypes[,135:188],na.rm=T),5:191]<-
  (2-genotypes[rowMeans(genotypes[,5:104],na.rm=T)>rowMeans(genotypes[,135:188],na.rm=T),5:191])

# Extract core region
genotypes<-genotypes[genotypes$POS>600000 & genotypes$POS<1000000,]

# Replace NA by -9
genotypes[is.na(genotypes)] = -9

# Function to replace up to 10 NA values by adjacent value if both adjacent values are equal
replNA<-function(x){
  
  # get runs
  conseq<-as.data.frame(cbind("length"=rle(x)$lengths,"value"=rle(x)$values,
                              "start"=c(1,cumsum(rle(x)$lengths)[-length(rle(x)$lengths)]+1)))
  conseq$end<-conseq$start+conseq$length-1
    
  naRuns<-conseq[conseq$value==-9,]
  naRuns<-naRuns[naRuns$start>1 & naRuns$end<length(x),] # remove first and last index if present
  
  # Check if the adjacent runs have the same value
  replace<-naRuns[which(x[naRuns$start-1]==x[naRuns$end+1]),]

  # get all positions to replace
  replacePos<-unlist(mapply(FUN = function(a, b) {
    seq(from = a, to = b, by = 1)
  }, a = replace$start, b = replace$end))
  
  # get values to replace them with
  replaceVal<-unlist(mapply(FUN = function(a, b) {
    rep(a, times = b)
  }, a = x[replace$start-1], b = replace$length))  
  
  # Replace these NA entries by their adjacent values
  x[replacePos]<-replaceVal
  
  return(x)
}

# Apply the NA replacement function multiple times
genotypes[,5:length(genotypes)]<-lapply(X = genotypes[,5:length(genotypes)],FUN = replNA)
genotypes[,5:length(genotypes)]<-lapply(X = genotypes[,5:length(genotypes)],FUN = replNA)
genotypes[,5:length(genotypes)]<-lapply(X = genotypes[,5:length(genotypes)],FUN = replNA)

# Function to replace up to 5 genotypes if they boarder to runs of the same value of at least 10 genotypes on either side
fixErrorsRegion<-function(x){

  # get runs
  conseq<-as.data.frame(cbind("length"=rle(x)$lengths,"value"=rle(x)$values,
                              "start"=c(1,cumsum(rle(x)$lengths)[-(length(rle(x)$lengths))]+1)))
  conseq$end<-conseq$start+conseq$length-1

  # Get short runs of up to 5 genotypes
  potErr<-conseq[conseq$length<6,]
  potErr<-potErr[potErr$start>1 & potErr$end<length(x),] # remove first and last index if present

  # Add the mean value of the region 10 bp north and south of the error region
  potErr$regionMean<-unlist(mapply(FUN = function(a, b) {
    indices<-(max(0,a-10)):min(b+10,max(conseq$end))  # get indices of region 10 bp north and south of the error region
    indices<-indices[!indices%in%c(a:b)] # remove indices of actual errors
    values<-x[indices] # get genotype values for these indices
    values[values==(-9)]<-NA # set -9 to NA
    if(length(which(is.na(values)))>15){meanVal=NA}else{meanVal<-mean(values,na.rm = T)}
    return(meanVal)
  }, a = potErr$start, b = potErr$end))

  # Remove potErr entries where the regionMean is NA
  potErr<-potErr[!is.na(potErr$regionMean),]

  # Remove regions that are unclear with what they would be replaced
  potErr<-potErr[potErr$regionMean<=0.1|(potErr$regionMean>=0.4 & potErr$regionMean<=0.6)|potErr$regionMean>=0.9,]

  # Get replacement value
  if(length(potErr$regionMean)>0) potErr$replaceVal<-plyr::round_any(potErr$regionMean,accuracy=1)

  # get all positions to replace
  replacePos<-unlist(mapply(FUN = function(a, b) {
    seq(from = a, to = b, by = 1)
  }, a = potErr$start, b = potErr$end))

  # get values to replace them with
  replaceVal<-unlist(mapply(FUN = function(a, b) {
    rep(a, times = b)
  }, a = potErr$replaceVal, b = potErr$length))

  # Replace these likely errors by their adjacent values
  x[replacePos]<-replaceVal

  return(x)
}

# Apply the error correction function multiple times
genotypes[,5:length(genotypes)]<-lapply(X=genotypes[,5:length(genotypes)],FUN=fixErrorsRegion)
genotypes[,5:length(genotypes)]<-lapply(X=genotypes[,5:length(genotypes)],FUN=fixErrorsRegion)
genotypes[,5:length(genotypes)]<-lapply(X=genotypes[,5:length(genotypes)],FUN=fixErrorsRegion)

# Move the recombinants to the end of the dataframe
recombinants<-c(80,48,2,17,99)
genotypes<-cbind(genotypes[,-(recombinants+4)],genotypes[,(recombinants+4)])

# Read in melpomene genes
melpGenes<-read.table("D:/Dropbox/Heliconius/HybridZones/geneAnnotations/Hmel2.5.gff3",sep="\t")
melpmRNA<-melpGenes[melpGenes$V3=="mRNA",]
optixMelpPos<-melpmRNA[grepl("HMEL001028-RA",melpmRNA$V9),]
melmRNAChr18<-melpmRNA[melpmRNA$V1=="Hmel218003o" & 
                         melpmRNA$V4>rangeMelOptix[1] & 
                         melpmRNA$V5<rangeMelOptix[2],]

# Other loci coordinates
otherOptix<-rbind(
cbind(772587,779428,"Wallbank_Ray"),
cbind(815837,816077,"Wallbank_Dennis_1"),
cbind(813673,814333,"Wallbank_Dennis_2"),
cbind(814176,817339,"JoeHanly_dennis"),
cbind(773272,787408,"JoeHanly_ray"),
cbind(718819,730976,"JoeHanly_band1"),
cbind(780304,796765,"JoeHanly_band2")
)

library(qdapTools)

# Function for making colours transparent for background shading
makeTransparent<-function(someColor, alpha=150)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

# function to get colour runs and plot them as background shading
plotGenotypes <- function(x, b,width=0.5) {
  colRuns<-as.data.frame(cbind("length"=as.numeric(rle(x)$lengths),
                               "value"=as.character(rle(x)$values),
                               "start"=c(1,cumsum(as.numeric(rle(x)$lengths))[-(length(rle(x)$lengths))]+1)))
  colRuns$startPos<-genotypes$POS[as.numeric(as.character(colRuns$start))]
  colRuns$endPos<-genotypes$POS[as.numeric(as.character(colRuns$start))+as.numeric(as.character(colRuns$length))-1]
  
  # Plot colour runs
  rect(xleft = colRuns$startPos,xright=colRuns$endPos,
       ybottom = b-width,ytop=b+width,
       col=makeTransparent(colRuns$value),border = NA)
  
  # Plot grey between colour runs
  if(length(colRuns$length)>1)
    rect(xleft = colRuns$endPos[-length(colRuns$endPos)]+1,xright=colRuns$startPos[-1]-1,
         ybottom = b-width,ytop=b+width,
         col=makeTransparent("grey"),border = NA)
}


# Convert genotypes to colours
colours<-t(as.matrix(genotypes[,5:191]))
colours[colours==0]<-"red"
colours[colours==1]<-"orange"
colours[colours==2]<-"yellow"
colours[colours==(-9)]<-"grey"

# Plot the optix haplotypes
par(mfrow=c(1,1),mar=c(4,4,1,1),cex.lab=1.8,cex.axis=1.5, xaxs="r",yaxs="r")
plot(max(genotypes$POS),length(genotypes),ylim=c(0,length(genotypes)+30),
     xlim=c(680000,850000),pch=NA,yaxt="n",
     xlab="positions [Mb]",ylab="",xaxt="n")
points(optix$win_mid,optix$LRT_mean/5+192,cex=0.5,type="l")
abline(h=192,lwd=3)
axis(1,at=seq(0,5e6,by=0.25e5),labels=seq(0,5,by=0.025))
rect(xleft = optixMelpPos$V4,
     xright=optixMelpPos$V5,
     ybottom=200,ytop=205,col="#C0202A",border = "#C0202A")
rect(xleft=as.integer(otherOptix[1:3,1]),
     xright=as.integer(otherOptix[1:3,2]),
     ybottom=208,ytop=212,col="blue")
rect(xleft=as.integer(otherOptix[4:6,1]),
     xright=as.integer(otherOptix[4:6,2]),
     ybottom=212,ytop=216,col="cornflowerblue")

indices<-c(1:95,97:127,129:181,183:185,187:191)
for(i in 1:length(colours[,1])){
  plotGenotypes(x=as.character(colours[i,]),b=indices[i])
}

# Convert colours to data frame
colours<-as.data.frame(colours)
names(colours)<-genotypes$POS[1:length(colours)]

# 4 dennis/ray no red in forewing band
points(rep(genotypes$POS[1:length(colours)],each=5),
       y=rep(187:191,times=length(colours[1,])),
       pch=15,col=as.vector(unlist(colours[183:187,])),cex=0.4)
abline(h=186,lwd=2)

# 3 dennis no ray
points(rep(genotypes$POS[1:length(colours)],each=3),
     y=rep(183:185,times=length(colours[1,])),
     pch=15,col=as.vector(unlist(colours[180:182,])),cex=0.4)
abline(h=182,lwd=2)

# plesseni-looking at optix
points(rep(genotypes$POS[1:length(colours)],each=53),
       y=rep(129:181,times=length(colours[1,])),
       pch=15,col=as.vector(unlist(colours[127:179,])),cex=0.4)
abline(h=128,lwd=2)

# phenotypically heterozygotes
points(rep(genotypes$POS[1:length(colours)],each=31),
       y=rep(97:127,times=length(colours[1,])),
       pch=15,col=as.vector(unlist(colours[96:126,])),cex=0.4)
abline(h=96,lwd=2)

# malleti
points(rep(genotypes$POS[1:length(colours)],each=95),
       y=rep(1:95,times=length(colours[1,])),
       pch=15,col=as.vector(unlist(colours[1:95,])),cex=0.4)


###########################################################################

# Plot haplotypes of the region including cortex and domeless/washout

###########################################################################


# Read in the vcf file of sites that have an allele frequency difference of at least 0.9 between plesseni and malleti-like individuals
genotypes<-read.table("cortex.diagnostic.SNPs.0.9hom.max5het.vcf",header=T,sep="\t")

# Specify plesseni allele as 1
genotypes[rowMeans(genotypes[,5:36],na.rm=T)<rowMeans(genotypes[,64:121],na.rm=T),5:length(genotypes)]<-
  (2-genotypes[rowMeans(genotypes[,5:36],na.rm=T)<rowMeans(genotypes[,64:121],na.rm=T),5:length(genotypes)])

# Replace NA by -9
genotypes[is.na(genotypes)] = -9

# Replace up to 10 NA values by adjacent value if both adjacent values are equal
genotypes[,5:length(genotypes)]<-lapply(X = genotypes[,5:length(genotypes)],FUN = replNA)

# Replace up to 5 genotypes if they boarder to runs of the same value
genotypes[,5:length(genotypes)]<-lapply(X=genotypes[,5:length(genotypes)],FUN=fixErrorsRegion)
genotypes[,5:length(genotypes)]<-lapply(X=genotypes[,5:length(genotypes)],FUN=fixErrorsRegion)
genotypes[,5:length(genotypes)]<-lapply(X=genotypes[,5:length(genotypes)],FUN=fixErrorsRegion)

# Convert genotypes to colours
colours<-t(as.matrix(genotypes[,5:length(genotypes)]))
colours[colours==0]<-"red"
colours[colours==1]<-"orange"
colours[colours==2]<-"yellow"
colours[colours==(-9)]<-"grey"

# Plot haplotypes
par(mfrow=c(1,1),mar=c(4,4,1,1),cex.lab=1.8,cex.axis=1.5, xaxs="r",yaxs="r")
plot(max(genotypes$POS),length(genotypes),ylim=c(0,length(genotypes)+5),
     pch=NA,yaxt="n",xlim=c(min(genotypes$POS),max(genotypes$POS)),
     xlab="positions [Mb]",ylab="",xaxt="n")
axis(1,at=seq(0,5e6,by=1e5),labels=seq(0,5,by=0.1))

# Plot genes cortex and dome/washout
rect(xleft = 1413776,
     xright=1533113,
     ybottom=-4,ytop=-1,col="grey",border = NA) # cortex
rect(xleft = 1684714 ,
     xright=1718943,
     ybottom=-4,ytop=-1,col="#C0202A",border = "#C0202A")

# X axis coordinates for haplotypes
indices<-c(1:32,34:60,62:119,122:122) 

# Plot the haplotypes
for(i in 1:length(colours[,1])){
  plotGenotypes(x=as.character(as.matrix(colours[i,])),b=indices[i])
}

# Convert colours to data frame
colours<-as.data.frame(colours)
names(colours)<-genotypes$POS[1:length(colours)]

# yellow spot (lowland-like), red distributed / short dennis bar (highland-like)
points(rep(genotypes$POS[1:length(colours)],each=1),
       y=rep(122,times=length(colours[1,])),
       pch=15,col=as.vector(unlist(colours[118,])),cex=0.4)
abline(h=120.5,lwd=2)
abline(h=123.5,lwd=2)

# highland-like in both traits
points(rep(genotypes$POS[1:length(colours)],each=58),
       y=rep(62:119,times=length(colours[1,])),
       pch=15,col=as.vector(unlist(colours[60:117,])),cex=0.4)
abline(h=61,lwd=2)

# no yellow spot (like highland) but red scales distributed like lowland
points(rep(genotypes$POS[1:length(colours)],each=27),
       y=rep(34:60,times=length(colours[1,])),
       pch=15,col=as.vector(unlist(colours[33:59,])),cex=0.4)
abline(h=33,lwd=2)

# lowland-like
points(rep(genotypes$POS[1:length(colours)],each=32),
       y=rep(1:32,times=length(colours[1,])),
       pch=15,col=as.vector(unlist(colours[1:32,])),cex=0.4)

