# R code for Figure 3

##################################################################################
# Read in input files
##################################################################################

library(stringr)

# Read in erato scaffold info to get scaffold length for stitching scaffolds together
ref.scaff <- read.table('Heliconius_erato_demophoon_v1_-_scaffolds.fa.fai')
ref.scaff<-ref.scaff[,1:2]
names(ref.scaff)<-c("scaff","length")
ref.scaff$add<-c(0,cumsum(ref.scaff$length)[-length(ref.scaff$scaff)])
ref.scaff$CHR<-as.integer(substr(ref.scaff$scaff,7,8))

# Get mid positions for chromosomes
chr<-aggregate(ref.scaff[,3],list(ref.scaff$CHR),min)
names(chr)<-c("CHR","add")
chr$length<-aggregate(ref.scaff[,2],list(ref.scaff$CHR),sum)[,2]
chr$mid<-chr$add+chr$length/2

# Read in Fst values of erato
erato<-read.table("erato.fst.10kb.idx.withHMM",header=T,sep="\t")
erato$win_mid_add<-erato$midPos+ref.scaff[match(as.character(erato$chr),as.character(ref.scaff$scaff)),"add"] # make a column with additive positions
erato<-erato[!is.na(erato$win_mid_add),]
erato$CHR = as.factor(as.integer(str_sub(erato$chr,7,8)))
erato<-erato[!is.na(erato$CHR),]

# Read in Fst values of melpomene with erato coordinates
melpomene<-read.table("fst.10kb.idx.withHMM.eratoLiftOver",header=T,sep="\t")
melpomene<-melpomene[order(melpomene$scaffold,as.integer(melpomene$pos)),]
melpomene$win_mid_add<-melpomene$pos+ref.scaff[match(as.character(melpomene$scaffold),as.character(ref.scaff$scaff)),"add"] # make a column with additive positions
melpomene<-melpomene[!is.na(melpomene$win_mid_add),]
melpomene$CHR = as.factor(as.integer(str_sub(melpomene$scaffold,7,8)))

# Read in the melpomene GWAS files with erato positions:

# number of forewing bands
wntAMelp<-read.table("gwas.wntA.p-val.pos.eratoLiftOver",header=T)
wntAMelp$win_mid_add<-wntAMelp$Position+ref.scaff[match(as.character(wntAMelp$scaffold),as.character(ref.scaff$scaff)),"add"]

# red patterns
optixMelp<-read.table("gwas.optix.NGSadmixK3.Asso2.cov.lrt0.gz.10kb.p-val.pos.eratoLiftOver",header=T)
optixMelp$win_mid_add<-optixMelp$Position+ref.scaff[match(as.character(optixMelp$scaffold),as.character(ref.scaff$scaff)),"add"]

# distribution of red scales
Nmelp<-read.table("gwas.N.NGSadmixK3.Asso2.model3.cov.lrt0.gz.10kb.p-val.pos.eratoLiftOver",header=T)
Nmelp$win_mid_add<-Nmelp$Position+ref.scaff[match(as.character(Nmelp$scaffold),as.character(ref.scaff$scaff)),"add"]

# yellow forewing base spot
ybsMelp<-read.table("gwas.ybs.NGSadmixK3.Asso2.model3.cov.lrt0.gz.10kb.p-val.pos.eratoLiftOver",header=T)
ybsMelp$win_mid_add<-ybsMelp$Position+ref.scaff[match(as.character(ybsMelp$scaffold),as.character(ref.scaff$scaff)),"add"]


# Read in the erato GWAS file containing all loci:

# Number of forewing bands
wntAErato<-read.table("gwas.wntA.cov.lrt0.gz.10kb.p-val.pos",header=T)
wntAErato$win_mid_add<-wntAErato$win_mid+ref.scaff[match(as.character(wntAErato$Chromosome),as.character(ref.scaff$scaff)),"add"]

# Red patterns
optixErato<-read.table("gwas.optix.cov.lrt0.gz.10kb.p-val.pos",header=T)
optixErato$win_mid_add<-optixErato$win_mid+ref.scaff[match(as.character(optixErato$Chromosome),as.character(ref.scaff$scaff)),"add"]

# Shape of forewing band edge
RoErato<-read.table("gwas.ro.model3.cov.lrt0.gz.10kb.p-val.pos",header=T)
RoErato$win_mid_add<-RoErato$win_mid+ref.scaff[match(as.character(RoErato$Chromosome),as.character(ref.scaff$scaff)),"add"]

# Yellow forewing base spot
ybsErato<-read.table("gwas.ybs.model3.cov.lrt0.gz.10kb.p-val.pos",header=T)
ybsErato$win_mid_add<-ybsErato$win_mid+ref.scaff[match(as.character(ybsErato$Chromosome),as.character(ref.scaff$scaff)),"add"]


# set extremely tiny p-values which gave infinite when -log10 tranformed to 250
wntAErato$logMinP[is.infinite(wntAErato$logMinP)]<-250
optixErato$logMinP[is.infinite(optixErato$logMinP)]<-250
RoErato$logMinP[is.infinite(RoErato$logMinP)]<-250
ybsErato$logMinP[is.infinite(ybsErato$logMinP)]<-250


############################################################################################
# Make plot for Figure 3A
############################################################################################


# Get erato FST peak positions:
er10<-erato[erato$chr=="Herato1001",]
wntAMid<-er10$win_mid_add[which.max(er10$Fst)]
er13<-erato[erato$chr=="Herato1301",]
nMid<-er13$win_mid_add[which.max(er13$Fst)]
er15<-erato[erato$chr=="Herato1505",]
roMid<-er15$win_mid_add[which.max(er15$Fst)]
er18<-erato[erato$chr=="Herato1801",]
optixMid<-er18$win_mid_add[which.max(er18$Fst)]

# Mid positions of erato peaks for yellow lines
midPos<-cbind(optixMid,wntAMid,roMid,nMid)


# Compute smoothed -log10 P-values for Fig 3A
{
  layout(matrix(c(1:5)),
         widths=c(1,1,1,1,1), heights=c(1,1,0.1,1,1))
  par(mar=c(0,0.5,0,0.5),oma=c(5,7,3,1),
      cex.lab=2.3,cex.axis=2,mgp=c(4,1,0),family='sans')
  
  # Plot Fst erato
  plot(erato$win_mid_add,erato$Fst,col=NA,xaxs="i", yaxs="i",
       xlab="",yaxt="n",xaxt="n",ylab="",ylim=c(0,1),xlim=c(0,382825554))
  axis(3,at=chr$mid[seq(1,21,by=2)],labels = c(chr$CHR[-21],"Z")[seq(1,21,by=2)],tick = T,cex=2.3,xpd=T)
  axis(3,at=chr$mid[seq(2,21,by=2)],labels = c(chr$CHR[-21],"Z")[seq(2,21,by=2)],tick = T,cex=2.3,xpd=T)
  rect(chr$add,0,chr$length+chr$add,1,
       col=rep(c("white","lightgrey"),n=21),border = NA)
  # abline(v=midPos,col="yellow",lwd=3)
  points(erato$win_mid_add,erato$Fst,col="grey50",cex=0.7,pch=16)
  points(erato$win_mid_add,erato$Fst,col=ifelse(erato$HMM>1,"black",NA),cex=0.7,pch=16)
  axis(2,at=c(0,0.5,1),labels = c(0.00,0.50,1.00),las=2)
  
  # Plot GWAS erato
  plot(wntAErato$win_mid_add,wntAErato$logMinP,col=NA,xaxs="i", yaxs="i",
       xlab="",yaxt="n",xaxt="n",ylab="",ylim=c(0,300),xlim=c(0,382825554))
  rect(chr$add,0,chr$length+chr$add,300,
       col=rep(c("white","lightgrey"),n=21),border = NA)
  # abline(v=midPos,col="yellow",lwd=3)
  smEra<-supsmu(wntAErato$win_mid_add,wntAErato$logMinP, span=0.0001)
  lines(smEra$x, smEra$y,col="#CCCC00",xaxs="i", yaxs="i",
        xlab="",xaxt="n",cex=2,pch=16)
  smEra<-supsmu(RoErato$win_mid_add,RoErato$logMinP, span=0.0001)
  lines(smEra$x, smEra$y,col="#0000CC",xaxs="i", yaxs="i",
         xlab="",xaxt="n",cex=2,pch=16)
  smEra<-supsmu(optixErato$win_mid_add,optixErato$logMinP, span=0.0001)
  lines(smEra$x, smEra$y,col="#C0202A",xaxs="i", yaxs="i",
         xlab="",xaxt="n",cex=2,pch=16)
  smEra<-supsmu(ybsErato$win_mid_add,ybsErato$logMinP, span=0.0001)
  lines(smEra$x, smEra$y,col="black",xaxs="i", yaxs="i",
         xlab="",xaxt="n",cex=2,pch=16)
  axis(2,at=seq(0,260,by=50),labels = seq(0,260,by=50),las=2)
  
  plot.new()
  
  # Plot Fst melpomene
  plot(melpomene$win_mid_add,melpomene$Fst,col=NA,xaxs="i", yaxs="i",
       xlab="",yaxt="n",xaxt="n",ylab="",ylim=c(0,1),xlim=c(0,382825554))
  rect(chr$add,0,chr$length+chr$add,1,
       col=rep(c("white","lightgrey"),n=21),border = NA)
  # abline(v=midPos,col="yellow",lwd=3)
  points(melpomene$win_mid_add,melpomene$Fst,col="grey50",cex=0.7,pch=16)
  points(melpomene$win_mid_add,melpomene$Fst,col=ifelse(melpomene$HMM>1,"black",NA),cex=0.7,pch=16)
  axis(2,at=c(0,0.5),labels = c(0.00,0.50),las=2)
  
  
  # Plot GWAS melpomene
  plot(wntAMelp$win_mid_add,wntAMelp$logMinP,col=NA,xaxs="i", yaxs="i",
       xlab="",yaxt="n",xaxt="n",ylab="",ylim=c(0,25),xlim=c(0,382825554))
  rect(chr$add,0,chr$length+chr$add,25,
       col=rep(c("white","lightgrey"),n=21),border = NA)
  # abline(v=midPos,col="yellow",lwd=3)
  smMel<-supsmu(wntAMelp$win_mid_add,wntAMelp$logMinP, span=0.001)
  lines(smMel$x,smMel$y,col="#CCCC00",xaxs="i", yaxs="i",
         xlab="",xaxt="n",cex=2,pch=16)
  smMel<-supsmu(Nmelp$win_mid_add,Nmelp$logMinP, span=0.001)
  lines(smMel$x,smMel$y,col="#CC33FF",xaxs="i", yaxs="i",
         xlab="",xaxt="n",cex=2,pch=16)
  smMel<-supsmu(optixMelp$win_mid_add,optixMelp$logMinP, span=0.001)
  lines(smMel$x,smMel$y,col="#C0202A",xaxs="i", yaxs="i",
         xlab="",xaxt="n",cex=2,pch=16)
  smMel<-supsmu(ybsMelp$win_mid_add,ybsMelp$logMinP, span=0.001)
  lines(smMel$x,smMel$y,col="black",xaxs="i", yaxs="i",
         xlab="",xaxt="n",cex=2,pch=16)
  axis(2,at=c(0,10,20),labels = c(0,10,20),las=2)
}


############################################################################################
# Make plot for Figure 3B: Zoomed in to the four major peak regions for H. erato
############################################################################################

# Read in pi values
thetaLativ<-read.table("PC062_merged.PL.AD.HAPCUT2.inf0.05.lativitta.10kb.windowed.pi",header=T)
thetaNota<-read.table("PC062_merged.PL.AD.HAPCUT2.inf0.05.notabilis.10kb.windowed.pi",header=T)
names(thetaLativ)<-c("Chr","start","end","n","pi")
names(thetaNota)<-c("Chr","start","end","n","pi")
thetaLativ$WinCenter<-thetaLativ$start+5000
thetaNota$WinCenter<-thetaNota$start+5000

thetaLativ50kb<-read.table("PC062_merged.PL.AD.HAPCUT2.inf0.05.lativitta.50kb.windowed.pi",header=T)
thetaNota50kb<-read.table("PC062_merged.PL.AD.HAPCUT2.inf0.05.notabilis.50kb.windowed.pi",header=T)
names(thetaLativ50kb)<-c("Chr","start","end","n","pi")
names(thetaNota50kb)<-c("Chr","start","end","n","pi")
thetaLativ50kb$WinCenter<-thetaLativ50kb$start+25000
thetaNota50kb$WinCenter<-thetaNota50kb$start+25000

# Extract wntA scaffold
thetaLativWntA<-thetaLativ[thetaLativ$Chr=="Herato1001",]
thetaNotaWntA<-thetaNota[thetaNota$Chr=="Herato1001",]
thetaLativ50kbWntA<-thetaLativ50kb[thetaLativ50kb$Chr=="Herato1001",]
thetaNota50kbWntA<-thetaNota50kb[thetaNota50kb$Chr=="Herato1001",]

# Extract Ro scaffold
thetaLativRo<-thetaLativ[thetaLativ$Chr=="Herato1301",]
thetaNotaRo<-thetaNota[thetaNota$Chr=="Herato1301",]
thetaLativRo50kb<-thetaLativ50kb[thetaLativ50kb$Chr=="Herato1301",]
thetaNotaRo50kb<-thetaNota50kb[thetaNota50kb$Chr=="Herato1301",]

# Extract Optix scaffold
thetaLativOptix<-thetaLativ[thetaLativ$Chr=="Herato1801",]
thetaNotaOptix<-thetaNota[thetaNota$Chr=="Herato1801",]
thetaLativ50kbOptix<-thetaLativ50kb[thetaLativ50kb$Chr=="Herato1801",]
thetaNota50kbOptix<-thetaNota50kb[thetaNota50kb$Chr=="Herato1801",]

# Extract Optix scaffold
thetaLativN<-thetaLativ[thetaLativ$Chr=="Herato1505",]
thetaNotaN<-thetaNota[thetaNota$Chr=="Herato1505",]
thetaLativ50kbN<-thetaLativ50kb[thetaLativ50kb$Chr=="Herato1505",]
thetaNota50kbN<-thetaNota50kb[thetaNota50kb$Chr=="Herato1505",]

# Read in gene annotations
eratoGenes<-read.table("heliconius_erato_demophoon_v1_core_32_85_1_sort_uniq.gff",sep="\t")
eratomRNA<-eratoGenes[eratoGenes$V3=="mRNA",]
optixEratoPos<-eratomRNA[grepl("evm.TU.Herato1801.64",eratomRNA$V9),]
cortexEratoPos<-eratomRNA[grepl("evm.model.Herato1505.85",eratomRNA$V9),]
wntAEratoPos<-eratomRNA[grepl("evm.model.Herato1001.160",eratomRNA$V9),]
vvlEratoPos<-eratomRNA[grepl("evm.model.Herato1301.536|evm.model.Herato1301.537",eratomRNA$V9),]
rsp3EratoPos<-eratomRNA[grepl("evm.model.Herato1301.538",eratomRNA$V9),]
vvlRsp3EratoPos<-rbind(vvlEratoPos,rsp3EratoPos)

# Get erato peak positions:
er10<-erato[erato$chr=="Herato1001",]
midWntA<-er10$midPos[which.max(er10$Fst)]
er13<-erato[erato$chr=="Herato1301",]
midRo<-er13$midPos[which.max(er13$Fst)]
er15<-erato[erato$chr=="Herato1505",]
midN<-er15$midPos[which.max(er15$Fst)]
er18<-erato[erato$chr=="Herato1801",]
midOptix<-er18$midPos[which.max(er18$Fst)]

# iHS 
iHSerato<-read.table("PC062_merged_all.erato.notabilis_lativitta.iHS.out",header=T)
omegaErato<-read.table("PC062_merged_all.omega.notabilis_lativitta.wHeader.out",header=T)


maxGWAS<-max(wntAErato$logMinP ,optixErato$logMinP ,RoErato$logMinP ,ybsErato$logMinP )+20

# Make Fig3B
{
  
  #Make a plot of all of them:
  layout(matrix(c(1:24), 6, 4, byrow = F),
         widths=c(1,1,1,1), heights=c(0.4,1,1,1))
  par(mar=c(0,0.5,0,0.5),oma=c(5,7,1,1),cex.lab=2.3,cex.axis=2,
      mgp=c(4,1,0))
  
  # wntA region
  # get the range of positions to plot
  rangeEraWntA<-c(midWntA-1e6,midWntA+1e6) 
  
  # get the genes for these positions
  ermRNAChr10<-eratomRNA[eratomRNA$V1=="Herato1001" & 
                           eratomRNA$V4>rangeEraWntA[1] & 
                           eratomRNA$V5<rangeEraWntA[2],]
  
  # Set up the plotting area
  plot(0,y=0.3,xlim=rangeEraWntA,xaxt="n",yaxt="n",pch=25,ylim=c(0,1),xlab="")

  # Plot the genes
  rect(xleft = ermRNAChr10$V4,
       xright=ermRNAChr10$V5,
       ybottom=0,ytop=0.3,col="grey",border = NA)
  rect(xleft = wntAEratoPos$V4,
       xright=wntAEratoPos$V5,
       ybottom=0,ytop=0.3,col="#C0202A",border = NA)
  points(mean(c(wntAEratoPos$V4,wntAEratoPos$V5)),y=0.3,pch=25,col="#C0202A",bg="#C0202A")
  text(x = wntAEratoPos$V4,y=0.7,labels = expression(italic("WntA")),cex=2)
  
  # Plot FST
  plot(erato$midPos[erato$chr=="Herato1001"],erato$Fst[erato$chr=="Herato1001"],
       cex=0.7,ylab="Fst",col=ifelse(erato$HMMstate>1,"#C0202A","grey50"),pch=NA,xlim=rangeEraWntA,xaxt="n",yaxt="n",ylim=c(0,1))
  points(erato$midPos[erato$chr=="Herato1001"],erato$Fst[erato$chr=="Herato1001"],
         cex=0.7,xlim=rangeEraWntA,xaxt="n",yaxt="n",ylim=c(0,1),pch=19,
         col=ifelse(erato$HMMstate[erato$chr=="Herato1001"]>1,"#C0202A","grey50"))
  axis(2,at=seq(0,1,by=0.5),labels=seq(0,1,by=0.5),las=2)
  
  # Plot pi
  plot(thetaLativWntA$WinCenter,thetaLativWntA$pi-thetaNotaWntA$pi,ylim=c(-0.02,0.02),
       xlim=rangeEraWntA,cex=0.7,ylab="",col=NA,yaxt="n",xaxt="n")
  axis(2,at=seq(-0.01,0.01,by=0.02),labels=seq(-0.01,0.01,by=0.02),las=2)
  points(thetaLativWntA$WinCenter,thetaNotaWntA$pi-thetaLativWntA$pi,
         cex=0.7,col="grey",pch=19)#,col=ifelse((thetaLativWntA$pi-thetaNotaWntA$pi)<0,"#C0202A","#CCCC00"))
  points(thetaLativ50kbWntA$WinCenter,thetaNota50kbWntA$pi-thetaLativ50kbWntA$pi,type="l",
         cex=0.7,col="black",pch=19,lwd=2)
  abline(h=0,lwd=0.5)
  
  # Plot iHS
  plot(iHSerato$START[iHSerato$CHROM=="Herato1001"]+5000,
       iHSerato$iHS_NOTABILIS[iHSerato$CHROM=="Herato1001"],
       xlim=rangeEraWntA,xaxt="n",yaxt="n",ylim=c(0,1),
       type="l",col="gold2",lwd=2)
  points(iHSerato$START[iHSerato$CHROM=="Herato1001"]+5000,
         iHSerato$iHS_LATIVITTA[iHSerato$CHROM=="Herato1001"],
         xlim=rangeEraWntA,type="l",col="red",lwd=2)
  axis(2,at=seq(0,0.5,by=0.5),labels=seq(0,0.5,by=0.5),las=2)
  
  # plot omega
  plot(omegaErato$POS_NOTABILIS[omegaErato$CHROM=="Herato1001"]+5000,
       omegaErato$Omega_MAX_NOTABILIS[omegaErato$CHROM=="Herato1001"],
       xlim=rangeEraWntA,xaxt="n",yaxt="n",lwd=2,ylim=c(0,5000),
       type="l",col="gold2")
  points(omegaErato$POS_NOTABILIS[omegaErato$CHROM=="Herato1001"]+5000,
         omegaErato$Omega_MAX_LATIVITTA[omegaErato$CHROM=="Herato1001"],
         xlim=rangeEraWntA,lwd=2,type="l",col="red")
  axis(2,at=seq(0,7000,by=3000),labels=seq(0,7000,by=3000),las=2)
  
  # Plot GWAS results
  plot(wntAErato$win_mid[wntAErato$Chromosome=="Herato1001"],wntAErato$logMinP[wntAErato$Chromosome=="Herato1001"],pch=NA,
       cex=0.7,ylab="LRT",xlim=rangeEraWntA,xaxt="n",yaxt="n",ylim=c(0,maxGWAS),col="black")
  points(optixErato$win_mid[optixErato$Chromosome=="Herato1001"],optixErato$logMinP[optixErato$Chromosome=="Herato1001"],
         cex=0.7,col="#C0202A",pch=19)
  points(RoErato$win_mid[RoErato$Chromosome=="Herato1001"],RoErato$logMinP[RoErato$Chromosome=="Herato1001"],
         cex=0.7,col="#0000CC",pch=19)
  points(ybsErato$win_mid[RoErato$Chromosome=="Herato1001"],ybsErato$logMinP[RoErato$Chromosome=="Herato1001"],
         cex=0.7,col="black",pch=19)
  points(wntAErato$win_mid[wntAErato$Chromosome=="Herato1001"],wntAErato$logMinP[wntAErato$Chromosome=="Herato1001"],
         cex=0.7,xlim=rangeEraWntA,xaxt="n",yaxt="n",ylim=c(0,480),col="#CCCC00",pch=19)
  axis(2,at=seq(0,500,by=200),labels=seq(0,500,by=200),las=2)
  axis(1,at=seq(0,5e8,by=5e5),labels = seq(0,5e2,by=0.5))
  axis(1,at=midWntA,labels = "Herato1001 (Mb)",outer=T,line=2,tick=F,cex.axis=2.2)
  
  # Ro region (same approach as described above for wntA)
  rangeEraRo<-c(midRo-1e6,midRo+1e6)
  ermRNAChr13<-eratomRNA[eratomRNA$V1=="Herato1301" & 
                           eratomRNA$V4>rangeEraRo[1] & 
                           eratomRNA$V5<rangeEraRo[2],]
  plot(0,y=0.3,xlim=rangeEraRo,xaxt="n",yaxt="n",pch=25,ylim=c(0,1),xlab="")
  # abline(v=midRo,col="yellow",lwd=3)
  rect(xleft = ermRNAChr13$V4,
       xright=ermRNAChr13$V5,
       ybottom=0,ytop=0.3,col="grey",border = NA)
  rect(xleft = vvlRsp3EratoPos$V4,
       xright=vvlRsp3EratoPos$V5,
       ybottom=0,ytop=0.3,col="#C0202A",border = NA)
  points(x=rowMeans(cbind(vvlRsp3EratoPos$V4,vvlRsp3EratoPos$V5)),
         y=c(0.3,0.3,0.3),pch=25,col="#C0202A",bg="#C0202A")
  text(x = vvlRsp3EratoPos$V4[2:3],y=c(0.65,0.55),pos=c(2,4),
       labels = c(expression(italic("vvl")),expression(italic("rsp3"))),cex=2)
  plot(erato$midPos[erato$chr=="Herato1301"],erato$Fst[erato$chr=="Herato1301"],
       cex=0.7,ylab="",xlim=rangeEraRo,xaxt="n",yaxt="n",ylim=c(0,1),pch=NA)
  #abline(v=midRo,col="yellow",lwd=3)
  points(erato$midPos[erato$chr=="Herato1301"],erato$Fst[erato$chr=="Herato1301"],
         cex=0.7,ylab="",xlim=rangeEraRo,xaxt="n",yaxt="n",ylim=c(0,1),pch=19,
         col=ifelse(erato$HMMstate[erato$chr=="Herato1301"]>1,"#C0202A","grey50"))
  plot(thetaLativRo$WinCenter,thetaNotaRo$pi-thetaLativRo$pi,ylim=c(-0.02,0.02),xaxt="n",
       xlim=rangeEraRo,cex=0.7,ylab="",col="grey",yaxt="n",pch=19,yaxt="n")
  #abline(v=midRo,col="yellow",lwd=3)
  points(thetaLativRo$WinCenter,thetaNotaRo$pi-thetaLativRo$pi,ylim=c(-0.02,0.02),xaxt="n",
         xlim=rangeEraRo,cex=0.7,ylab="",col="grey",yaxt="n",pch=19,yaxt="n")
  points(thetaLativRo50kb$WinCenter,thetaNotaRo50kb$pi-thetaLativRo50kb$pi,ylim=c(-0.02,0.02),xaxt="n",
         xlim=rangeEraRo,cex=0.7,ylab="",col="black",yaxt="n",pch=19,yaxt="n",type="l",lwd=2)
  abline(h=0,lwd=0.5)
  plot(iHSerato$START[iHSerato$CHROM=="Herato1301"]+5000,
       iHSerato$iHS_NOTABILIS[iHSerato$CHROM=="Herato1301"],
       xlim=rangeEraRo,xaxt="n",yaxt="n",ylim=c(0,1),
       type="l",col="gold2",lwd=2)
  points(iHSerato$START[iHSerato$CHROM=="Herato1301"]+5000,
         iHSerato$iHS_LATIVITTA[iHSerato$CHROM=="Herato1301"],
         xlim=rangeEraRo,type="l",col="red",lwd=2)
  plot(omegaErato$POS_NOTABILIS[omegaErato$CHROM=="Herato1301"]+5000,
       omegaErato$Omega_MAX_NOTABILIS[omegaErato$CHROM=="Herato1301"],
       xlim=rangeEraRo,xaxt="n",yaxt="n",lwd=2,ylim=c(0,5000),
       type="l",col="gold2")
  points(omegaErato$POS_NOTABILIS[omegaErato$CHROM=="Herato1301"]+5000,
         omegaErato$Omega_MAX_LATIVITTA[omegaErato$CHROM=="Herato1301"],
         xlim=rangeEraRo,lwd=2,type="l",col="red")
  plot(RoErato$win_mid[RoErato$Chromosome=="Herato1301"],RoErato$logMinP[RoErato$Chromosome=="Herato1301"],
       cex=0.7,ylab="",xlim=rangeEraRo,xaxt="n",yaxt="n",ylim=c(0,maxGWAS))
  #abline(v=midRo,col="yellow",lwd=3)
  points(optixErato$win_mid[optixErato$Chromosome=="Herato1301"],optixErato$logMinP[optixErato$Chromosome=="Herato1301"],
         cex=0.7,col="#C0202A",pch=19)
  points(wntAErato$win_mid[wntAErato$Chromosome=="Herato1301"],wntAErato$logMinP[wntAErato$Chromosome=="Herato1301"],
         cex=0.7,xlim=rangeEraRo,xaxt="n",yaxt="n",ylim=c(0,480),col="#CCCC00",pch=19)
  points(ybsErato$win_mid[RoErato$Chromosome=="Herato1301"],ybsErato$logMinP[RoErato$Chromosome=="Herato1301"],
         cex=0.7,col="black",pch=19)
  points(RoErato$win_mid[RoErato$Chromosome=="Herato1301"],RoErato$logMinP[RoErato$Chromosome=="Herato1301"],
         cex=0.7,col="#0000CC",pch=19)
  axis(1,at=seq(0,5e8,by=1e6),labels = seq(0,5e2,by=1))
  axis(1,at=seq(5e5,5e8,by=1e6),labels = seq(0.5,5e2,by=1))
  axis(1,at=midRo,labels = "Herato1301 (Mb)",outer=T,line=2,tick=F,cex.axis=2.2)
  
  # Cortex region (same approach as described above for wntA)
  rangeEraN<-c(midN-1e6,midN+1e6)
  ermRNAChr15<-eratomRNA[eratomRNA$V1=="Herato1505" & 
                           eratomRNA$V4>rangeEraN[1] & 
                           eratomRNA$V5<rangeEraN[2],]
  plot(0,y=0.3,xlim=rangeEraN,xaxt="n",yaxt="n",pch=NA,ylim=c(0,1),xlab="")
  # abline(v=midN,col="yellow",lwd=3)
  rect(xleft = ermRNAChr15$V4,
       xright=ermRNAChr15$V5,
       ybottom=0,ytop=0.3,col="grey",border = NA)
  rect(xleft = cortexEratoPos$V4,
       xright=cortexEratoPos$V5,
       ybottom=0,ytop=0.3,col="#C0202A",border = NA)
  points(mean(c(cortexEratoPos$V4,cortexEratoPos$V5)),y=0.3,pch=25,col="#C0202A",bg="#C0202A")
  text(x = cortexEratoPos$V4,y=0.7,labels = expression(italic("cortex")),cex=2)
  plot(erato$midPos[erato$chr=="Herato1505"],erato$Fst[erato$chr=="Herato1505"],
       cex=0.7,ylab="Fst",xlim=rangeEraN,xaxt="n",yaxt="n",ylim=c(0,1),pch=NA)
  #abline(v=midN,col="yellow",lwd=3)
  points(erato$midPos[erato$chr=="Herato1505"],erato$Fst[erato$chr=="Herato1505"],
         cex=0.7,ylab="",xlim=rangeEraN,xaxt="n",yaxt="n",ylim=c(0,1),pch=19,
         col=ifelse(erato$HMMstate[erato$chr=="Herato1505"]>1,"#C0202A","grey50"))
  plot(thetaLativN$WinCenter,thetaLativN$pi-thetaNotaN$pi,ylim=c(-0.02,0.02),xaxt="n",
       xlim=rangeEraN,cex=0.7,ylab="",col=NA,yaxt="n",pch=19)
  #abline(v=midN,col="yellow",lwd=3)
  points(thetaLativN$WinCenter,thetaNotaN$pi-thetaLativN$pi,
         cex=0.7,col="grey",pch=19)
  points(thetaLativ50kbN$WinCenter,thetaNota50kbN$pi-thetaLativ50kbN$pi,type="l",
         cex=0.7,col="black",pch=19,lwd=2)
  abline(h=0,lwd=0.5)
  plot(iHSerato$START[iHSerato$CHROM=="Herato1505"]+5000,
       iHSerato$iHS_NOTABILIS[iHSerato$CHROM=="Herato1505"],
       xlim=rangeEraN,xaxt="n",yaxt="n",ylim=c(0,1),
       type="l",col="gold2",lwd=2)
  points(iHSerato$START[iHSerato$CHROM=="Herato1505"]+5000,
         iHSerato$iHS_LATIVITTA[iHSerato$CHROM=="Herato1505"],
         xlim=rangeEraN,type="l",col="red",lwd=2)
  plot(omegaErato$POS_NOTABILIS[omegaErato$CHROM=="Herato1505"]+5000,
       omegaErato$Omega_MAX_NOTABILIS[omegaErato$CHROM=="Herato1505"],
       xlim=rangeEraN,xaxt="n",yaxt="n",lwd=2,ylim=c(0,5000),
       type="l",col="gold2")
  points(omegaErato$POS_NOTABILIS[omegaErato$CHROM=="Herato1505"]+5000,
         omegaErato$Omega_MAX_LATIVITTA[omegaErato$CHROM=="Herato1505"],
         xlim=rangeEraN,lwd=2,type="l",col="red")
  plot(RoErato$win_mid[RoErato$Chromosome=="Herato1505"],RoErato$logMinP[RoErato$Chromosome=="Herato1505"],
       cex=0.7,ylab="",xlim=rangeEraN,xaxt="n",yaxt="n",col="black",pch=NA,ylim=c(0,maxGWAS))
  #abline(v=midN,col="yellow",lwd=3)
  points(optixErato$win_mid[optixErato$Chromosome=="Herato1505"],optixErato$logMinP[optixErato$Chromosome=="Herato1505"],
         cex=0.7,col="#C0202A",pch=19)
  points(wntAErato$win_mid[wntAErato$Chromosome=="Herato1505"],wntAErato$logMinP[wntAErato$Chromosome=="Herato1505"],
         cex=0.7,col="#CCCC00",pch=19)
  points(RoErato$win_mid[RoErato$Chromosome=="Herato1505"],RoErato$logMinP[RoErato$Chromosome=="Herato1505"],
         cex=0.7,col="#0000CC",pch=19)
  points(ybsErato$win_mid[RoErato$Chromosome=="Herato1505"],ybsErato$logMinP[RoErato$Chromosome=="Herato1505"],
         cex=0.7,col="black",pch=19)
  axis(1,at=seq(0,5e8,by=5e5),labels = seq(0,5e2,by=0.5))
  axis(1,at=midN,labels = "Herato1505 (Mb)",outer=T,line=2,tick=F,cex.axis=2.2)
  
  
  # Optix region (same approach as described above for wntA)
  rangeEraOptix<-c(midOptix-1e6,midOptix+1e6)
  ermRNAChr18<-eratomRNA[eratomRNA$V1=="Herato1801" & 
                           eratomRNA$V4>rangeEraOptix[1] & 
                           eratomRNA$V5<rangeEraOptix[2],]
  plot(0,y=0.3,xlim=rangeEraOptix,xaxt="n",yaxt="n",pch=NA,ylim=c(0,1),xlab="")
  #abline(v=midOptix,col="yellow",lwd=3)
  rect(xleft = ermRNAChr18$V4,
       xright=ermRNAChr18$V5,
       ybottom=0,ytop=0.3,col="grey",border = NA)
  rect(xleft = optixEratoPos$V4,
       xright=optixEratoPos$V5,
       ybottom=0,ytop=0.3,col="#C0202A",border = NA)
  points(mean(c(optixEratoPos$V4,optixEratoPos$V5)),y=0.3,pch=25,col="#C0202A",bg="#C0202A")
  text(x = optixEratoPos$V4,y=0.7,labels = expression(italic("optix")),cex=2)
  plot(erato$midPos[erato$chr=="Herato1801"],erato$Fst[erato$chr=="Herato1801"],
       cex=0.7,ylab="Fst",xlim=rangeEraOptix,xaxt="n",yaxt="n",ylim=c(0,1),pch=NA)
  #abline(v=midOptix,col="yellow",lwd=3)
  points(erato$midPos[erato$chr=="Herato1801"],erato$Fst[erato$chr=="Herato1801"],
         cex=0.7,ylab="",xlim=rangeEraOptix,xaxt="n",yaxt="n",ylim=c(0,1),pch=19,
         col=ifelse(erato$HMMstate[erato$chr=="Herato1801"]>1,"#C0202A","grey50"))
  plot(thetaLativOptix$WinCenter,thetaLativOptix$pi-thetaNotaOptix$pi,ylim=c(-0.02,0.02),xaxt="n",
       xlim=rangeEraOptix,cex=0.7,ylab="",col=NA,yaxt="n")
  # abline(v=midOptix,col="yellow",lwd=3)
  points(thetaLativOptix$WinCenter,thetaNotaOptix$pi-thetaLativOptix$pi,
         cex=0.7,col="grey",pch=19)
  points(thetaLativ50kbOptix$WinCenter,thetaNota50kbOptix$pi-thetaLativ50kbOptix$pi,type="l",
         cex=0.7,col="black",pch=19,lwd=2)
  abline(h=0,lwd=0.5)
  plot(iHSerato$START[iHSerato$CHROM=="Herato1801"]+5000,
       iHSerato$iHS_NOTABILIS[iHSerato$CHROM=="Herato1801"],
       xlim=rangeEraOptix,xaxt="n",yaxt="n",ylim=c(0,1),
       type="l",col="gold2",lwd=2)
  points(iHSerato$START[iHSerato$CHROM=="Herato1801"]+5000,
         iHSerato$iHS_LATIVITTA[iHSerato$CHROM=="Herato1801"],
         xlim=rangeEraOptix,type="l",col="red",lwd=2)
  plot(omegaErato$POS_NOTABILIS[omegaErato$CHROM=="Herato1801"]+5000,
       omegaErato$Omega_MAX_NOTABILIS[omegaErato$CHROM=="Herato1801"],
       xlim=rangeEraOptix,xaxt="n",yaxt="n",lwd=2,ylim=c(0,5000),
       type="l",col="gold2")
  points(omegaErato$POS_NOTABILIS[omegaErato$CHROM=="Herato1801"]+5000,
         omegaErato$Omega_MAX_LATIVITTA[omegaErato$CHROM=="Herato1801"],
         xlim=rangeEraOptix,lwd=2,type="l",col="red")
  plot(optixErato$win_mid[optixErato$Chromosome=="Herato1801"],optixErato$logMinP[optixErato$Chromosome=="Herato1801"],
       cex=0.7,ylab="",xlim=rangeEraOptix,xaxt="n",yaxt="n",pch=NA,ylim=c(0,maxGWAS))
  #abline(v=midOptix,col="yellow",lwd=3)
  points(wntAErato$win_mid[wntAErato$Chromosome=="Herato1801"],wntAErato$logMinP[wntAErato$Chromosome=="Herato1801"],
         cex=0.7,col="#CCCC00",pch=19)
  points(RoErato$win_mid[RoErato$Chromosome=="Herato1801"],RoErato$logMinP[RoErato$Chromosome=="Herato1801"],
         cex=0.7,col="#0000CC",pch=19)
  points(optixErato$win_mid[optixErato$Chromosome=="Herato1801"],optixErato$logMinP[optixErato$Chromosome=="Herato1801"],
         cex=0.7,col="#C0202A",pch=19)
  axis(1,at=seq(0,5e8,by=5e5),labels = seq(0,5e2,by=0.5))
  axis(1,at=midOptix,labels = "Herato1801 (Mb)",outer=T,line=2,tick=F,cex.axis=2.2) 
}


############################################################################################
# Make plot for Figure 3C: Zoomed in to the four major peak regions for H. melpomene
############################################################################################

# Read in pi values
thetaMalle<-read.table("D:/Dropbox/Heliconius/HybridZones/melpomene/Theta/run164_merged.PL.AD.HAPCUT2.inf0.05.malleti.10kb.windowed.pi",header=T)
thetaPless<-read.table("D:/Dropbox/Heliconius/HybridZones/melpomene/Theta/run164_merged.PL.AD.HAPCUT2.inf0.05.plesseni.10kb.windowed.pi",header=T)
names(thetaMalle)<-c("Chr","start","end","n","pi")
names(thetaPless)<-c("Chr","start","end","n","pi")
thetaMalle$WinCenter<-thetaMalle$start+5000
thetaPless$WinCenter<-thetaPless$start+5000

thetaMalle50kb<-read.table("D:/Dropbox/Heliconius/HybridZones/melpomene/Theta/run164_merged.PL.AD.HAPCUT2.inf0.05.malleti.50kb.windowed.pi",header=T)
thetaPless50kb<-read.table("D:/Dropbox/Heliconius/HybridZones/melpomene/Theta/run164_merged.PL.AD.HAPCUT2.inf0.05.plesseni.50kb.windowed.pi",header=T)
names(thetaMalle50kb)<-c("Chr","start","end","n","pi")
names(thetaPless50kb)<-c("Chr","start","end","n","pi")
thetaMalle50kb$WinCenter<-thetaMalle50kb$start+25000
thetaPless50kb$WinCenter<-thetaPless50kb$start+25000

# Extract wntA scaffold
thetaMalleWntA<-thetaMalle[thetaMalle$Chr=="Hmel210001o",]
thetaPlessWntA<-thetaPless[thetaPless$Chr=="Hmel210001o",]
thetaMalle50kbWntA<-thetaMalle50kb[thetaMalle50kb$Chr=="Hmel210001o",]
thetaPless50kbWntA<-thetaPless50kb[thetaPless50kb$Chr=="Hmel210001o",]

# Extract Ro scaffold
thetaMalleRo<-thetaMalle[thetaMalle$Chr=="Hmel213001o",]
thetaPlessRo<-thetaPless[thetaPless$Chr=="Hmel213001o",]
thetaMalle50kbRo<-thetaMalle50kb[thetaMalle50kb$Chr=="Hmel213001o",]
thetaPless50kbRo<-thetaPless50kb[thetaPless50kb$Chr=="Hmel213001o",]

# Extract Optix scaffold
thetaMalleOptix<-thetaMalle[thetaMalle$Chr=="Hmel218003o",]
thetaPlessOptix<-thetaPless[thetaPless$Chr=="Hmel218003o",]
thetaMalle50kbOptix<-thetaMalle50kb[thetaMalle50kb$Chr=="Hmel218003o",]
thetaPless50kbOptix<-thetaPless50kb[thetaPless50kb$Chr=="Hmel218003o",]

# Extract N scaffold
thetaMalleN<-thetaMalle[thetaMalle$Chr=="Hmel215003o",]
thetaPlessN<-thetaPless[thetaPless$Chr=="Hmel215003o",]
thetaMalle50kbN<-thetaMalle50kb[thetaMalle50kb$Chr=="Hmel215003o",]
thetaPless50kbN<-thetaPless50kb[thetaPless50kb$Chr=="Hmel215003o",]

# GWAS
N<-read.table("D:gwas.N.NGSadmixK3.Asso2.model3.cov.lrt0.gz.10kb.p-val.pos",header=T)
N<-N[!is.na(N$Chromosome),]

ybs<-read.table("gwas.ybs.NGSadmixK3.Asso2.model3.cov.lrt0.gz.10kb.p-val.pos",header=T)
ybs<-ybs[!is.na(ybs$Chromosome),]

# optix quantitative trait
optix<-read.table("gwas.optix.NGSadmixK3.Asso2.cov.lrt0.gz.10kb.p-val.pos",header=T)
optix<-optix[!is.na(optix$Chromosome),]

# wntA quantitative trait
wntA<-read.table("gwas.wntA.NGSadmixK3.Asso2.cov.lrt0.gz.10kb.p-val.pos",header=T)
wntA<-wntA[!is.na(wntA$Chromosome),]

iHSmel<-read.table("run164_merged_all.melpomene.plesseni_malleti.iHS.out",header=T)
omegaMel<-read.table("run164_merged_all.omega.plesseni_malleti.wHeader.out",header=T)

library(stringr)

# Read in scaffold info:
ref.scaff <- read.table('Hmel2.5.fa.fai')
ref.scaff<-ref.scaff[,1:2]
names(ref.scaff)<-c("scaff","length")
ref.scaff$add<-c(0,cumsum(ref.scaff$length)[-length(ref.scaff$scaff)])
ref.scaff<-ref.scaff[ref.scaff$scaff!="Hmel_complete_mtDNA",] # remove mitochondrial DNA sites
ref.scaff$CHR<-as.integer(substr(ref.scaff$scaff,6,7))
ref.scaff<-ref.scaff[!grepl(ref.scaff$scaff,pattern = "Hmel200"),] # remove scaffolds not assigned to a chromosome

# Generate chromosome info:
chr<-aggregate(ref.scaff[,3],list(ref.scaff$CHR),min)
names(chr)<-c("CHR","add")
chr$length<-aggregate(ref.scaff[,2],list(ref.scaff$CHR),sum)[,2]
chr$mid<-chr$add+chr$length/2

# Read in Fst values
fst<-read.table("melpomene.fst.10kb.idx.withHMM",header=T,sep="\t")
fst<-fst[!grepl(fst$chr,pattern = "Hmel200"),]
fst$win_mid_add<-fst$midPos+ref.scaff[match(as.character(fst$chr),as.character(ref.scaff$scaff)),"add"]
fst<-fst[!is.na(fst$win_mid_add),]
fst$CHR = as.factor(as.integer(str_sub(fst$chr,6,7)))
melpomene<-fst
maxGWAS<-max(wntA$logMinP,optix$logMinP,N$logMinP,ybs$logMinP)
maxGWAS<-40

# Read in gene annotations
melpGenes<-read.table("Hmel2.5.gff3",sep="\t")
melpmRNA<-melpGenes[melpGenes$V3=="mRNA",]
optixMelpPos<-melpmRNA[grepl("HMEL001028-RA",melpmRNA$V9),]
cortexMelpPos<-melpmRNA[grepl("HMEL000025-RH",melpmRNA$V9),]
wntaMelpPos<-melpmRNA[grepl("HMEL018100g1.t1",melpmRNA$V9),]
vvlRsp3MelpPos<-melpmRNA[grepl("HMEL011784g1.t1|HMEL020781g1.t1",melpmRNA$V9),]


# Plot melpomene zoomed into the four regions (window min p-value) for panel C
{
  layout(matrix(c(1:24), 6, 4, byrow = F),
         widths=c(1,1,1,1), heights=c(0.4,1,1,1))
  par(mar=c(0,0.5,0,0.5),oma=c(5,7,1,1),cex.lab=2.3,cex.axis=2,mgp=c(4,1,0))
  
  # wntA
  mid<-mean(melpomene$midPos[melpomene$chr=="Hmel210001o"&melpomene$Fst>0.3])
  midWntA<-3335500
  rangeMelWntA<-c(midWntA-1e6,midWntA+1e6)
  melmRNAChr10<-melpmRNA[melpmRNA$V1=="Hmel210001o" & 
                           melpmRNA$V4>rangeMelWntA[1] & 
                           melpmRNA$V5<rangeMelWntA[2],]
  plot(0,y=0.3,xlim=rangeMelWntA,xaxt="n",yaxt="n",pch=25,ylim=c(0,1),xlab="")
  # abline(v=midWntA,col="yellow",lwd=3)
  rect(xleft = melmRNAChr10$V4,
       xright=melmRNAChr10$V5,
       ybottom=0,ytop=0.3,col="grey",border = NA)
  rect(xleft = wntaMelpPos$V4,
       xright=wntaMelpPos$V5,
       ybottom=0,ytop=0.3,col="#C0202A",border = NA)
  points(mean(c(wntaMelpPos$V4,wntaMelpPos$V5)),y=0.3,pch=25,col="#C0202A",bg="#C0202A")
  text(x = wntaMelpPos$V4,y=0.7,labels = expression(italic("WntA")),cex=2)
  plot(melpomene$midPos[melpomene$chr=="Hmel210001o"],melpomene$Fst[melpomene$chr=="Hmel210001o"],
       cex=0.7,ylab="",xlim=rangeMelWntA,xaxt="n",yaxt="n",pch=NA,ylim=c(0,1))
  # abline(v=midWntA,col="yellow",lwd=3)
  points(melpomene$midPos[melpomene$chr=="Hmel210001o"],melpomene$Fst[melpomene$chr=="Hmel210001o"],
         cex=0.7,ylab="Fst",xlim=rangeMelWntA,xaxt="n",yaxt="n",ylim=c(0,1),pch=19,
         col=ifelse(melpomene$HMMstate[melpomene$chr=="Hmel210001o"]>1,"#C0202A","grey50"))
  axis(2,at=seq(0,1,by=0.5),labels=seq(0,1,by=0.5),las=2)
  plot(thetaMalleWntA$WinCenter,thetaPlessWntA$pi-thetaMalleWntA$pi,ylim=c(-0.02,0.02),
       xlim=rangeMelWntA,cex=0.7,ylab=expression(paste(Delta, "pi pless-mall")),col=NA,yaxt="n",
       xlab="",xpd=T,xaxt="n",pch=19)
  # abline(v=midWntA,col="yellow",lwd=3)
  axis(2,at=seq(-0.01,0.01,by=0.02),labels=seq(-0.01,0.01,by=0.02),las=2)
  points(thetaMalleWntA$WinCenter,thetaPlessWntA$pi-thetaMalleWntA$pi,
         cex=0.7,col="grey",pch=19)
  points(thetaMalle50kbWntA$WinCenter,thetaPless50kbWntA$pi-thetaMalle50kbWntA$pi,type="l",
         cex=0.7,col="black",pch=19,lwd=2)
  abline(h=0,lwd=0.5)
  # abline(v=midWntA,col="yellow",lwd=3)
  plot(iHSmel$START[iHSmel$CHROM=="Hmel210001o"]+5000,
       iHSmel$iHS_PLESSENI[iHSmel$CHROM=="Hmel210001o"],
       xlim=rangeMelWntA,xaxt="n",yaxt="n",ylim=c(0,1),
       type="l",col="gold2",lwd=2)
  points(iHSmel$START[iHSmel$CHROM=="Hmel210001o"]+5000,
         iHSmel$iHS_MALLETI[iHSmel$CHROM=="Hmel210001o"],
         xlim=rangeMelWntA,type="l",col="red",lwd=2)
  axis(2,at=c(0,0.5),labels=c(0,0.5),las=2)
  plot(omegaMel$POS_PLESSENI[omegaMel$CHROM=="Hmel210001o"]+5000,
       omegaMel$Omega_MAX_PLESSENI[omegaMel$CHROM=="Hmel210001o"],
       xlim=rangeMelWntA,col="gold2",type="l",lwd=2,ylim=c(0,3800),
       xaxt="n",yaxt="n")
  points(omegaMel$POS_PLESSENI[omegaMel$CHROM=="Hmel210001o"]+5000,
         omegaMel$Omega_MAX_MALLETI[omegaMel$CHROM=="Hmel210001o"],
         xlim=rangeMelWntA,col="red",type="l",lwd=2)
  axis(2,at=c(0,2000),labels=c(0,2000),las=2)
  plot(wntA$win_mid,wntA$logMinP,pch=NA,
       cex=0.7,ylab="LRT",xlim=rangeMelWntA,xaxt="n",yaxt="n",ylim=c(0,maxGWAS),col="black")
  points(optix$win_mid[optix$Chromosome=="Hmel210001o"],optix$logMinP[optix$Chromosome=="Hmel210001o"],
         cex=0.7,ylab="",xlim=rangeMelWntA,xaxt="n",yaxt="n",col="#C0202A",pch=19)
  points(ybs$win_mid[N$Chromosome=="Hmel210001o"],ybs$logMinP[N$Chromosome=="Hmel210001o"],
         cex=0.7,ylab="",xlim=rangeMelWntA,xaxt="n",yaxt="n",col="black",pch=19)
  points(N$win_mid[N$Chromosome=="Hmel210001o"],N$logMinP[N$Chromosome=="Hmel210001o"],
         cex=0.7,ylab="",xlim=rangeMelWntA,xaxt="n",yaxt="n",col="#CC33FF",pch=19)
  points(wntA$win_mid[wntA$Chromosome=="Hmel210001o"],wntA$logMinP[wntA$Chromosome=="Hmel210001o"],
         cex=0.7,ylab="",xlim=rangeMelWntA,xaxt="n",yaxt="n",col="#CCCC00",pch=19)
  axis(2,at=seq(0,30,by=10),labels=seq(0,30,by=10),las=2)
  axis(1,at=midWntA,labels = "Hmel210001o (Mb)",outer=T,line=2,tick=F,cex.axis=2.2)
  axis(1,at=seq(0,5e8,by=5e5),labels = seq(0,5e2,by=0.5))
  
  
  # Ro
  mid<-mean(melpomene$midPos[melpomene$chr=="Hmel213001o"&melpomene$Fst>0.3])
  midRo<-10352500
  rangeMelRo<-c(midRo-1e6,midRo+1e6)
  melmRNAChr13<-melpmRNA[melpmRNA$V1=="Hmel213001o" & 
                           melpmRNA$V4>rangeMelRo[1] & 
                           melpmRNA$V5<rangeMelRo[2],]
  plot(0,y=0.3,xlim=rangeMelRo,xaxt="n",yaxt="n",pch=25,ylim=c(0,1),xlab="")
  # abline(v=midRo,col="yellow",lwd=3)
  rect(xleft = melmRNAChr13$V4,
       xright=melmRNAChr13$V5,
       ybottom=0,ytop=0.3,col="grey",border = NA)
  rect(xleft = vvlRsp3MelpPos$V4,
       xright=vvlRsp3MelpPos$V5,
       ybottom=0,ytop=0.3,col="#C0202A",border = NA)
  points(x=rowMeans(cbind(vvlRsp3MelpPos$V4,vvlRsp3MelpPos$V5)),
         y=c(0.3,0.3),pch=25,col="#C0202A",bg="#C0202A")
  text(x = vvlRsp3MelpPos$V4,y=c(0.65,0.55),pos=c(2,4),
       labels = c(expression(italic("vvl")),expression(italic("rsp3"))),cex=2)
  plot(melpomene$midPos[melpomene$chr=="Hmel213001o"],melpomene$Fst[melpomene$chr=="Hmel213001o"],
       cex=0.7,ylab="",xlim=rangeMelRo,xaxt="n",yaxt="n",pch=NA,ylim=c(0,1))
  # abline(v=midRo,col="yellow",lwd=3)
  points(melpomene$midPos[melpomene$chr=="Hmel213001o"],melpomene$Fst[melpomene$chr=="Hmel213001o"],
         cex=0.7,ylab="",xlim=rangeMelRo,xaxt="n",yaxt="n",ylim=c(0,1),pch=19,
         col=ifelse(melpomene$HMMstate[melpomene$chr=="Hmel213001o"]>1,"#C0202A","grey50"))
  plot(thetaMalleRo$WinCenter,thetaPlessRo$pi-thetaMalleRo$pi,ylim=c(-0.02,0.02),xaxt="n",yaxt="n",
       xlim=rangeMelRo,cex=0.7,ylab="",col="grey",yaxt="n",pch=19)
  # abline(v=midRo,col="yellow",lwd=3)
  points(thetaMalleRo$WinCenter,thetaPlessRo$pi-thetaMalleRo$pi,
         cex=0.7,col="grey",pch=19)
  points(thetaMalle50kbRo$WinCenter,thetaPless50kbRo$pi-thetaMalle50kbRo$pi,
         type="l",cex=0.7,col="black",lwd=2)
  abline(h=0,lwd=0.5)
  plot(iHSmel$START[iHSmel$CHROM=="Hmel213001o"]+5000,
       iHSmel$iHS_PLESSENI[iHSmel$CHROM=="Hmel213001o"],
       xlim=rangeMelRo,xaxt="n",yaxt="n",ylim=c(0,1),
       type="l",col="gold2",lwd=2)
  points(iHSmel$START[iHSmel$CHROM=="Hmel213001o"]+5000,
         iHSmel$iHS_MALLETI[iHSmel$CHROM=="Hmel213001o"],
         xlim=rangeMelRo,type="l",col="red",lwd=2)
  plot(omegaMel$POS_PLESSENI[omegaMel$CHROM=="Hmel213001o"]+5000,
       omegaMel$Omega_MAX_PLESSENI[omegaMel$CHROM=="Hmel213001o"],
       xlim=rangeMelRo,col="gold2",type="l",lwd=2,ylim=c(0,3800),
       xaxt="n",yaxt="n")
  points(omegaMel$POS_PLESSENI[omegaMel$CHROM=="Hmel213001o"]+5000,
         omegaMel$Omega_MAX_MALLETI[omegaMel$CHROM=="Hmel213001o"],
         xlim=rangeMelRo,col="red",type="l",lwd=2)
  plot(ro$win_mid,ro$logMinP,pch=NA,yaxt="n",
       cex=0.7,ylab="",xlim=rangeMelRo,xaxt="n",yaxt="n",ylim=c(0,maxGWAS))
  # abline(v=midRo,col="yellow",lwd=3)
  points(optix$win_mid[optix$Chromosome=="Hmel213001o"],optix$logMinP[optix$Chromosome=="Hmel213001o"],
         cex=0.7,ylab="",xlim=rangeMelRo,xaxt="n",yaxt="n",ylim=c(0,maxGWAS),col="#C0202A",pch=19)
  points(N$win_mid[N$Chromosome=="Hmel213001o"],N$logMinP[N$Chromosome=="Hmel213001o"],
         cex=0.7,ylab="",xlim=rangeMelRo,xaxt="n",yaxt="n",ylim=c(0,maxGWAS),col="#CC33FF",pch=19)
  points(ybs$win_mid[N$Chromosome=="Hmel213001o"],ybs$logMinP[N$Chromosome=="Hmel213001o"],
         cex=0.7,ylab="",xlim=rangeMelRo,xaxt="n",yaxt="n",ylim=c(0,maxGWAS),col="black",pch=19)
  points(wntA$win_mid[wntA$Chromosome=="Hmel213001o"],wntA$logMinP[wntA$Chromosome=="Hmel213001o"],
         cex=0.7,ylab="",xlim=rangeMelRo,xaxt="n",yaxt="n",ylim=c(0,maxGWAS),col="#CCCC00",pch=19)
  axis(1,at=midRo,labels = "Hmel213001o (Mb)",outer=T,line=2,tick=F,cex.axis=2.2)
  axis(1,at=seq(0.5e6,5e8,by=1e6),labels = seq(0.5,5e2,by=1))
  axis(1,at=seq(0,5e8,by=1e6),labels = seq(0,5e2,by=1))
  
  # cortex
  mid<-mean(melpomene$midPos[melpomene$chr=="Hmel215003o"&melpomene$Fst>0.3])
  midN<-1442500
  rangeMelN<-c(midN-1e6,midN+1e6)
  melmRNAChr15<-melpmRNA[melpmRNA$V1=="Hmel215003o" & 
                           melpmRNA$V4>rangeMelN[1] & 
                           melpmRNA$V5<rangeMelN[2],]
  plot(0,y=0.3,xlim=rangeMelN,xaxt="n",yaxt="n",pch=25,ylim=c(0,1),xlab="")
  # abline(v=midN,col="yellow",lwd=3)
  rect(xleft = melmRNAChr15$V4,
       xright=melmRNAChr15$V5,
       ybottom=0,ytop=0.3,col="grey",border = NA)
  rect(xleft = cortexMelpPos$V4,
       xright=cortexMelpPos$V5,
       ybottom=0,ytop=0.3,col="#C0202A",border = NA)
  text(x = cortexMelpPos$V4,y=0.7,labels = expression(italic("cortex")),cex=2)
  plot(melpomene$midPos[melpomene$chr=="Hmel215003o"],melpomene$Fst[melpomene$chr=="Hmel215003o"],
       cex=0.7,ylab="",xlim=rangeMelN,xaxt="n",yaxt="n",pch=NA,ylim=c(0,1))
  # abline(v=midN,col="yellow",lwd=3)
  points(melpomene$midPos[melpomene$chr=="Hmel215003o"],melpomene$Fst[melpomene$chr=="Hmel215003o"],
         cex=0.7,ylab="",xlim=rangeMelN,xaxt="n",yaxt="n",ylim=c(0,1),pch=19,
         col=ifelse(melpomene$HMMstate[melpomene$chr=="Hmel215003o"]>1,"#C0202A","grey50"))
  plot(thetaMalleN$WinCenter,thetaPlessN$pi-thetaMalleN$pi,ylim=c(-0.02,0.02),xaxt="n",
       xlim=rangeMelN,cex=0.7,ylab="",col=NA,yaxt="n")
  # abline(v=midN,col="yellow",lwd=3)
  points(thetaMalleN$WinCenter,thetaPlessN$pi-thetaMalleN$pi,
         cex=0.7,col="grey",pch=19)#,col=ifelse((thetaMalleOptix$pi-thetaPlessOptix$pi)<0,"#C0202A","#CCCC00"))
  points(thetaMalle50kbN$WinCenter,thetaPless50kbN$pi-thetaMalle50kbN$pi,type="l",
         cex=0.7,col="black",pch=19,lwd=2)
  abline(h=0,lwd=0.5)
  plot(iHSmel$START[iHSmel$CHROM=="Hmel215003o"]+5000,
       iHSmel$iHS_PLESSENI[iHSmel$CHROM=="Hmel215003o"],
       xlim=rangeMelN,xaxt="n",yaxt="n",ylim=c(0,1),
       type="l",col="gold2",lwd=2)
  points(iHSmel$START[iHSmel$CHROM=="Hmel215003o"]+5000,
         iHSmel$iHS_MALLETI[iHSmel$CHROM=="Hmel215003o"],
         xlim=rangeMelN,type="l",col="red",lwd=2)
  plot(omegaMel$POS_PLESSENI[omegaMel$CHROM=="Hmel215003o"]+5000,
       omegaMel$Omega_MAX_PLESSENI[omegaMel$CHROM=="Hmel215003o"],
       xlim=rangeMelN,col="gold2",type="l",lwd=2,ylim=c(0,3800),
       xaxt="n",yaxt="n")
  points(omegaMel$POS_PLESSENI[omegaMel$CHROM=="Hmel215003o"]+5000,
         omegaMel$Omega_MAX_MALLETI[omegaMel$CHROM=="Hmel215003o"],
         xlim=rangeMelN,col="red",type="l",lwd=2)
  plot(N$win_mid,N$logMinP,ylim=c(0,maxGWAS),
       cex=0.7,ylab="",xlim=rangeMelN,xaxt="n",yaxt="n",pch=NA)
  # abline(v=midN,col="yellow",lwd=3)
  points(optix$win_mid[optix$Chromosome=="Hmel215003o"],optix$logMinP[optix$Chromosome=="Hmel215003o"],
         cex=0.7,ylab="",xlim=rangeMelN,xaxt="n",yaxt="n",ylim=c(0,maxGWAS),col="#C0202A",pch=19)
  points(wntA$win_mid[wntA$Chromosome=="Hmel215003o"],wntA$logMinP[wntA$Chromosome=="Hmel215003o"],
         cex=0.7,ylab="",xlim=rangeMelN,xaxt="n",yaxt="n",ylim=c(0,maxGWAS),col="#CCCC00",pch=19)
  points(N$win_mid[N$Chromosome=="Hmel215003o"],N$logMinP[N$Chromosome=="Hmel215003o"],
         cex=0.7,ylab="",xlim=rangeMelN,xaxt="n",yaxt="n",ylim=c(0,maxGWAS),col="#CC33FF",pch=19)
  points(ybs$win_mid[ybs$Chromosome=="Hmel215003o"],ybs$logMinP[ybs$Chromosome=="Hmel215003o"],
         cex=0.7,ylab="",xlim=rangeMelN,xaxt="n",yaxt="n",ylim=c(0,maxGWAS),col="black",pch=9)
  axis(1,at=midN,labels = "Hmel215003o (Mb)",outer=T,line=2,tick=F,cex.axis=2.2)
  axis(1,at=seq(0,5e8,by=5e5),labels = seq(0,5e2,by=0.5))
  
  
  # Optix
  mid<-mean(melpomene$midPos[melpomene$chr=="Hmel218003o"&melpomene$Fst>0.3])
  midOptix<-787500
  midOptix<-822500
  rangeMelOptix<-c(midOptix-1e6,midOptix+1e6)
  melmRNAChr18<-melpmRNA[melpmRNA$V1=="Hmel218003o" & 
                           melpmRNA$V4>rangeMelOptix[1] & 
                           melpmRNA$V5<rangeMelOptix[2],]
  plot(0,y=0.3,xlim=rangeMelOptix,xaxt="n",yaxt="n",pch=NA,ylim=c(0,1),xlab="")
  # abline(v=midOptix,col="yellow",lwd=3)
  rect(xleft = melmRNAChr18$V4,
       xright=melmRNAChr18$V5,
       ybottom=0,ytop=0.3,col="grey",border = NA)
  rect(xleft = optixMelpPos$V4,
       xright=optixMelpPos$V5,
       ybottom=0,ytop=0.3,col="#C0202A",border = NA)
  points(mean(c(optixMelpPos$V4,optixMelpPos$V5)),y=0.3,pch=25,col="#C0202A",bg="#C0202A")
  text(x = optixMelpPos$V4,y=0.7,labels = expression(italic("optix")),cex=2)
  plot(melpomene$midPos[melpomene$chr=="Hmel218003o"],melpomene$Fst[melpomene$chr=="Hmel218003o"],
       cex=0.7,ylab="",xlim=rangeMelOptix,xaxt="n",yaxt="n",ylim=c(0,1),pch=NA)
  points(melpomene$midPos[melpomene$chr=="Hmel218003o"],melpomene$Fst[melpomene$chr=="Hmel218003o"],
         cex=0.7,ylab="",xlim=rangeMelOptix,xaxt="n",yaxt="n",ylim=c(0,1),pch=19,
         col=ifelse(melpomene$HMMstate[melpomene$chr=="Hmel218003o"]>1,"#C0202A","grey50"))
  plot(thetaMalleOptix$WinCenter,thetaPlessOptix$pi-thetaMalleOptix$pi,xaxt="n",ylim=c(-0.02,0.02),
       xlim=rangeMelOptix,cex=0.7,ylab="",col=NA,yaxt="n",pch=19)
  # abline(v=midOptix,col="yellow",lwd=3)
  points(thetaMalleOptix$WinCenter,thetaPlessOptix$pi-thetaMalleOptix$pi,
         cex=0.7,col="grey",pch=19)#,col=ifelse((thetaMalleOptix$pi-thetaPlessOptix$pi)<0,"#C0202A","#CCCC00"))
  points(thetaMalle50kbOptix$WinCenter,thetaPless50kbOptix$pi-thetaMalle50kbOptix$pi,type="l",
         cex=0.7,col="black",pch=19,lwd=2)
  abline(h=0,lwd=0.5)
  plot(iHSmel$START[iHSmel$CHROM=="Hmel218003o"]+5000,
       iHSmel$iHS_PLESSENI[iHSmel$CHROM=="Hmel218003o"],
       xlim=rangeMelOptix,xaxt="n",yaxt="n",ylim=c(0,1),
       type="l",col="gold2",lwd=2)
  points(iHSmel$START[iHSmel$CHROM=="Hmel218003o"]+5000,
         iHSmel$iHS_MALLETI[iHSmel$CHROM=="Hmel218003o"],
         xlim=rangeMelOptix,type="l",col="red",lwd=2)
  plot(omegaMel$POS_PLESSENI[omegaMel$CHROM=="Hmel218003o"]+5000,
       omegaMel$Omega_MAX_PLESSENI[omegaMel$CHROM=="Hmel218003o"],
       xlim=rangeMelOptix,col="gold2",type="l",lwd=2,ylim=c(0,3800),
       xaxt="n",yaxt="n")
  points(omegaMel$POS_PLESSENI[omegaMel$CHROM=="Hmel218003o"]+5000,
         omegaMel$Omega_MAX_MALLETI[omegaMel$CHROM=="Hmel218003o"],
         xlim=rangeMelOptix,col="red",type="l",lwd=2)
  plot(optix$win_mid,optix$logMinP,ylim=c(0,maxGWAS),pch=NA,
       cex=0.7,ylab="",xlim=rangeMelOptix,xaxt="n",yaxt="n")
  # abline(v=midOptix,col="yellow",lwd=3)
  points(N$win_mid[N$Chromosome=="Hmel218003o"],N$logMinP[N$Chromosome=="Hmel218003o"],
         cex=0.7,ylab="",xlim=rangeMelOptix,xaxt="n",yaxt="n",col="#CC33FF",pch=19)
  points(ybs$win_mid[N$Chromosome=="Hmel218003o"],ybs$logMinP[N$Chromosome=="Hmel218003o"],
         cex=0.7,ylab="",xlim=rangeMelOptix,xaxt="n",yaxt="n",col="#CC33FF",pch=19)
  points(wntA$win_mid[wntA$Chromosome=="Hmel218003o"],wntA$logMinP[wntA$Chromosome=="Hmel218003o"],
         cex=0.7,ylab="",xlim=rangeMelOptix,xaxt="n",yaxt="n",col="#CCCC00",pch=19)
  points(optix$win_mid[optix$Chromosome=="Hmel218003o"],optix$logMinP[optix$Chromosome=="Hmel218003o"],
         cex=0.7,ylab="",xlim=rangeMelOptix,xaxt="n",yaxt="n",col="#C0202A",pch=19)
  axis(1,at=midOptix,labels = "Hmel218003o (Mb)",outer=T,line=2,tick=F,cex.axis=2.2)
  axis(1,at=seq(0,5e8,by=5e5),labels = seq(0,5e2,by=0.5))
}



