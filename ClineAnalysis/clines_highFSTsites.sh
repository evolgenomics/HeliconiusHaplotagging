# Bash and HZAR R code to get the clines for the most highly differentiated sites of each 100 kbp window

#####################################################################
# Get the positions with the highest FST value per 100 kbp window
#####################################################################

# Use realSFS to output per site alpha and beta
realSFS fst print melpomene.fst.idx > melpomene.fst.tmp

# Add header and Fst column (alpha/beta)
awk 'BEGIN{print "Chrom\tPos\talpha\tbeta\tFst"}{if($3!=0) print $0"\t"$3/$4; else print $0"\t"0}' melpomene.fst.tmp > melpomene.fst.pos

# In R, get the highest value per 100 kbp window
R
fst<-read.table("melpomene.fst.pos",header=T)

maxFst<-data.frame()
for(scaff in levels(fst$Chrom)){
  scaffFst<-fst[fst$Chrom==scaff,]
  if(length(scaffFst$Pos>0)){
    for(i in seq(1,scaffFst$Pos[length(scaffFst$Pos)],by=100000)){
     wind<-scaffFst[scaffFst$Pos>i & scaffFst$Pos<(i+99999),]
     if(length(wind$Pos)>0)
     maxFst<-rbind(maxFst,wind[which.max(wind$Fst),])
    }
  }
}

# Exclude windows that do not have a single value above 0.4 
fst0.25<-maxFst[maxFst$Fst>0.4,c(1,2)]

# Write out the tables
write.table(fst0.25,file="melpomene.Fst0.4.pos",quote=F,row.names=F,sep="\t")
quit() # quit R

# Extract these positions from the vcf files
for i in *.vcf.gz
do 
  outfile=`basename $i`
  outfile=${outfile%.vcf.gz}.highFstPos
  vcftools --gzvcf $i --positions melpomene.Fst0.4.pos --recode --out $outfile
done

# Concatenate the extracts of all scaffolds:
bcftools concat -Oz -o fst100kbMax0.4_melpomene.vcf.gz *.vcf

for i in group_{1..26}.frq
do
 cut -f 5 $i | sed "s/{FREQ}/${i%.frq}/" > ${i%.frq}.tmp
done
paste *tmp > fst100kbMax0.4_melpomene.frq

# Add Fst
paste melpomene.Fst0.4.pos fst100kbMax0.4_melpomene.frq > melpomene.100kbMax.Fst0.4.frq



#####################################################################
# Compute the clines for these positions with HZAR (R code)
#####################################################################

R

# Read in the locus number from the command line (this was run as an Rscript that takes a number as a single argument)
args<-commandArgs(TRUE)
locus<-args[1]

# Read in the input files
freq<-read.table("melpomene.100kbMax.Fst0.4.frq",header=T)
samples<-read.table("SAMPLE_INFO_Melpomene.grouping.definitive.list",sep="\t", header=TRUE)

# Order the groups alphanumerically
library(gtools)
freq<-freq[,c(names(freq)[!grepl("group",names(freq))],
  mixedsort(names(freq[grepl(pattern="group",names(freq))])))]

# Make the clines from high to low
freq[rowMeans(cbind(freq$group_1,freq$group_2))<rowMeans(cbind(freq$group_16,freq$group_17)),grepl(names(freq),pattern="group")]<-
  (1-freq[rowMeans(cbind(freq$group_1,freq$group_2))<rowMeans(cbind(freq$group_16,freq$group_17)),grepl(names(freq),pattern="group")])

print(names(freq))

# Get the transect positions and samples sizes for all sites
transect<-c(6.02,11.82,18.97,28.34,29.25,29.41,32,32.48,36.39,37.14,43.59,44.89,48.07,51.02,56.82,59.22,60.38,67.24,73.85)
sampleSize<-c(15,16,12,5,8,7,7,8,15,6,15,7,6,5,15,5,13,6,16)


# allele frequency at a position 'locus'
p<-as.double(freq[locus,paste0("group_",1:19)])

print(p)

# Set the chain length. This value is the default setting in the package.
chainLength=1e5

## Make each model run off a separate seed
mainSeed=list(A=c(596,528,124,978,544,99),
  B=c(528,124,978,544,99,596),
  C=c(124,978,544,99,596,528))

# Load packages
require(doMC)
require(hzar)

# Prepare an object holding everything
melpomene <- list()
## Space to hold the observed data
melpomene$obs <- list();
## Space to hold the models to fit
melpomene$models <- list();
## Space to hold the compiled fit requests
melpomene$fitRs <- list();
## Space to hold the output data chains
melpomene$runs <- list();
## Space to hold the analysed data
melpomene$analysis <- list();


# Generate a hzar.obsData object
melpomene$obs <- hzar.doMolecularData1DPops(distance = transect, pObs = p, nEff = sampleSize)

## Make a helper function to define scaling, tails and name of models
melpomene.loadAdaAmodel <- function(scaling,tails,id=paste(scaling,tails,sep="."))
melpomene$models[[id]] <<- hzar.makeCline1DFreq(melpomene$obs, scaling, tails)

# Add six models
melpomene.loadAdaAmodel("fixed","none","modelI")
melpomene.loadAdaAmodel("fixed","both","modelII")
melpomene.loadAdaAmodel("fixed","mirror","modelIII")
melpomene.loadAdaAmodel("free" ,"none","modelIV")
melpomene.loadAdaAmodel("free" ,"both","modelV")
melpomene.loadAdaAmodel("free","mirror","modelVI")

## Modify all models to focus on the region where the observed
## data were collected. Restrain the MCMC optimization search range:
## Observations were between 6.020 and 81.225 km.
melpomene$models <- sapply(melpomene$models,
hzar.model.addBoxReq,0,85,simplify=FALSE)

## Compile each of the models to prepare for fitting
melpomene$fitRs$init <- sapply(melpomene$models,
  hzar.first.fitRequest.old.ML,obsData=melpomene$obs,
  verbose=FALSE,simplify=FALSE)

## Run each model for an initial chain
melpomene$runs$init <- list()
melpomene$runs$init$modelI <- hzar.doFit(melpomene$fitRs$init$modelI)
melpomene$runs$init$modelII <- hzar.doFit(melpomene$fitRs$init$modelII)
melpomene$runs$init$modelIII <- hzar.doFit(melpomene$fitRs$init$modelIII)
melpomene$runs$init$modelIV <- hzar.doFit(melpomene$fitRs$init$modelIV)
melpomene$runs$init$modelV <- hzar.doFit(melpomene$fitRs$init$modelV)
melpomene$runs$init$modelVI <- hzar.doFit(melpomene$fitRs$init$modelVI)

## Compile a new set of fit requests using the initial chains
melpomene$fitRs$chains <- lapply(melpomene$runs$init,hzar.next.fitRequest)

## Replicate each fit request 3 times, keeping the original
## seeds while switching to a new seed channel.
melpomene$fitRs$chains <- hzar.multiFitRequest(melpomene$fitRs$chains,
  each=3,baseSeed=NULL)

## Just to be thorough, randomize the initial value for each fit
for(i in 1:18){
 melpomene$fitRs$chains[[i]]$modelParam$init["center"]<-runif(1,0,85)
 melpomene$fitRs$chains[[i]]$modelParam$init["width"]<-runif(1,0,85)
 if(i>13) melpomene$fitRs$chains[[i]]$modelParam$init["pMin"]<-runif(1,0,1)
 if(i>13) melpomene$fitRs$chains[[i]]$modelParam$init["pMax"]<-runif(1,0,1)
 if(i%in%c(4:6,13:15)) melpomene$fitRs$chains[[i]]$modelParam$init["deltaL"]=runif(1,0,85)
 if(i%in%c(4:6,13:15)) melpomene$fitRs$chains[[i]]$modelParam$init["tauL"]=runif(1,0,1)
 if(i%in%c(4:6,13:15)) melpomene$fitRs$chains[[i]]$modelParam$init["deltaR"]=runif(1,0,1)
 if(i%in%c(4:6,13:15)) melpomene$fitRs$chains[[i]]$modelParam$init["tauR"]=runif(1,0,1)
 if(i%in%c(7:9,16:18)) melpomene$fitRs$chains[[i]]$modelParam$init["tauM"]=runif(1,0,1)
 if(i%in%c(7:9,16:18)) melpomene$fitRs$chains[[i]]$modelParam$init["deltaM"]=runif(1,0,1)
}

# Run a chain of 3 runs for every fit request
melpomene$runs$chains <- hzar.doChain.multi(melpomene$fitRs$chains,
  doPar=TRUE,inOrder=FALSE,count=3)

## Did modelI converge?
summary(do.call(mcmc.list,
  lapply(melpomene$runs$chains[1:3],function(x) hzar.mcmc.bindLL(x[[3]]))))

## Create a model data group for the null model (expected allele
## frequency independent of distance along cline) to include in analysis.
melpomene$analysis$initDGs <- list(nullModel = hzar.dataGroup.null(melpomene$obs))

## Create a model data group (hzar.dataGroup object) for each
## model from the initial runs.
melpomene$analysis$initDGs$modelI <-
  hzar.dataGroup.add(melpomene$runs$init$modelI)
melpomene$analysis$initDGs$modelII <-
  hzar.dataGroup.add(melpomene$runs$init$modelII)
melpomene$analysis$initDGs$modelIII <-
  hzar.dataGroup.add(melpomene$runs$init$modelIII)
melpomene$analysis$initDGs$modelIV <-
  hzar.dataGroup.add(melpomene$runs$init$modelIV)
melpomene$analysis$initDGs$modelV <-
  hzar.dataGroup.add(melpomene$runs$init$modelV)
melpomene$analysis$initDGs$modelVI <-
  hzar.dataGroup.add(melpomene$runs$init$modelVI)

## Create a hzar.obsDataGroup object from the hzar.dataGroups just created
melpomene$analysis$oDG <-hzar.make.obsDataGroup(melpomene$analysis$initDGs)

# Use the same model names
melpomene$analysis$oDG <- hzar.copyModelLabels(melpomene$analysis$initDGs,
  melpomene$analysis$oDG)

## Convert all runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object.
melpomene$analysis$oDG <-hzar.make.obsDataGroup(lapply(melpomene$runs$chains,
  hzar.dataGroup.add),melpomene$analysis$oDG)

## Do model selection based on the AICc scores
print(melpomene$analysis$AICcTable <-
  hzar.AICc.hzar.obsDataGroup(melpomene$analysis$oDG));

## Print out the model with the minimum AICc score
print(melpomene$analysis$model.name <-
  rownames(melpomene$analysis$AICcTable)[[which.min(melpomene$analysis$AICcTable$AICc)]])

# Compute delta AICc (should be above 2)
deltaAICc<-sort(melpomene$analysis$AICcTable$AICc)[2]-sort(melpomene$analysis$AICcTable$AICc)[1]

## Extract the hzar.dataGroup object for the selected model
melpomene$analysis$model.selected <-
  melpomene$analysis$oDG$data.groups[[melpomene$analysis$model.name]]

## Print the maximum likelihood cline for the selected model
#print(hzar.get.ML.cline(melpomene$analysis$model.selected))

# Skip the following steps if the null model is best
if(melpomene$analysis$model.name!="nullModel"){

# Get the confidence intervals for the width and center:
confI<-hzar.qScores.dataGroup(melpomene$analysis$model.selected,probs = c(0, 0.25, 0.5, 0.75, 1))

# Get the range of values within 2 log-likelihoods of the ML value
# 2-unit support envelope allows visualization of uncertainty in model fit (Devitt et al. 2011; Macholan et al. 2011).
SI2<-hzar.getLLCutParam(melpomene$analysis$model.selected, c("width","center"), cutValue = 2)

# Get the locus name
locName<-paste0(as.character(freq[locus,"Chrom"]),"_",as.character(freq[locus,"Pos"]))

# Plot the cline
png(width=900, height=900, res=200, family="Arial",
  filename=paste0(locName,".png"),pointsize=8)

# plot the fitted cline
plot(transect,p,type="n",pch=21,bg=NA,xlab="Transect (km)", ylab="Allele freq.",ylim=c(0,1),
     main=paste(as.character(freq[locus,"Chrom"]),freq[locus,"Pos"],melpomene$analysis$model.name,
                 "c:",round(as.double(melpomene$analysis$model.selected$ML.cline$param.all)[1],2),
                 "w:",round(as.double(melpomene$analysis$model.selected$ML.cline$param.all)[2],2)))
hzar.plot.fzCline(melpomene$analysis$model.selected,add=TRUE,fzCol=rgb(100,100,100,30,maxColorValue=255), lwd=1,pch=NA)
hzar.plot.cline(melpomene$analysis$model.selected,add=TRUE,col="grey50", lty=2, lwd=2,cex=1,pch=16)
lines(rep(as.double(melpomene$analysis$model.selected$ML.cline$param.free)[1],2),
      c(mean(as.double(melpomene$analysis$model.selected$ML.cline$param.all)[3:4])-0.05,
      mean(as.double(melpomene$analysis$model.selected$ML.cline$param.all)[3:4])+0.05), lwd=2,col="grey50")

# wntA centre
abline(v=47.1 ,col="blue")

# Optix centre
abline(v=31.9,col="red")

dev.off()


# Write out the results
cat(c("AIC",locName,deltaAICc,melpomene$analysis$AICcTable$AICc,"\n"),
  file=paste0(locName,".txt"))

cat(c("bestModel",locName,melpomene$analysis$model.name,
 names(melpomene$analysis$model.selected$ML.cline$param.all),
 as.double(melpomene$analysis$model.selected$ML.cline$param.all),"\n"),
  file=paste0(locName,".txt"),append=T)

cat(c("CI",locName,confI$width,confI$center,"\n"),
  file=paste0(locName,".txt"),append=T)

cat(c("SI2",locName,as.double(SI2),"\n"),
  file=paste0(locName,".txt"),append=T)

# If the null model is best:
}else{
# Plot the cline
png(width=900, height=900, res=200, family="Arial",
  filename=paste0(locName,".png"),pointsize=8)


# plot the fitted cline
plot(transect,p,type="n",pch=21,bg=NA,xlab="Transect (km)", ylab="Allele freq.",ylim=c(0,1),
     main=paste(as.character(freq[locus,"Chrom"]),freq[locus,"Pos"],melpomene$analysis$model.name,
                 "c:",round(as.double(melpomene$analysis$model.selected$ML.cline$param.all)[1],2),
                 "w:",round(as.double(melpomene$analysis$model.selected$ML.cline$param.all)[2],2)))
hzar.plot.fzCline(melpomene$analysis$model.selected,add=TRUE,fzCol=rgb(100,100,100,30,maxColorValue=255), lwd=1,pch=NA)
hzar.plot.cline(melpomene$analysis$model.selected,add=TRUE,col="grey50", lty=2, lwd=2,cex=1,pch=16)
lines(rep(as.double(melpomene$analysis$model.selected$ML.cline$param.free)[1],2),
      c(mean(as.double(melpomene$analysis$model.selected$ML.cline$param.all)[3:4])-0.05,
      mean(as.double(melpomene$analysis$model.selected$ML.cline$param.all)[3:4])+0.05), lwd=2,col="grey50")

# wntA centre
abline(v=47.1 ,col="blue")

# Optix centre
abline(v=31.9,col="red")

dev.off()

# Print some info (skip all the confidence interval and the like)
cat(c("AIC",locName,deltaAICc,melpomene$analysis$AICcTable$AICc,"\n"),
  file=paste0(locName,".txt"))
cat(c("bestModel",locName,melpomene$analysis$model.name),file=paste0(locName,".txt"),append=T)
}



#####################################################################
# Compute the clines for these positions with HZAR (R code)
#####################################################################


# melpomene

# Check out the distribution of clines form modeling:
# Read in the clines
clinesMel<-read.table("D:/Dropbox/Heliconius/HybridZones/melpomene/clines/modelI_modelIV.list",header=T)
clinesMelIII<-read.table("D:/Dropbox/Heliconius/HybridZones/melpomene/clines/modelIII.list",header=T)
transect<-c(6.02,11.82,18.97,28.34,29.25,29.41,32,32.48,36.39,37.14,43.59,44.89,48.07,51.02,56.82,59.22,60.38,67.24,73.85)
sampleSize<-c(15,16,12,5,8,7,7,8,15,6,15,7,6,5,15,5,13,6,16)


# plot the fitted clinesMel of all high Fst values
par(mfrow=c(1,1))
plot(transect,rep(1,times=length(transect)),type="n",pch=21,bg=NA,xlab="Transect (km)", ylab="Allele frequency", ylim=c(0,1),
     main="clinesMel of highest Fst site per 100kb window")

for(i in 1:length(clinesMel$Chrom)){
  curve(clinesMel$pMax[i] + (clinesMel$pMin[i] - clinesMel$pMax[i]) * (1/(1 + exp(-((x - clinesMel$center[i]) * 4/clinesMel$width[i])))),
        add=T,col=rgb(0.5,0.5,0.5,alpha = 0.5,maxColorValue = 1))
}

for(i in 1:length(clinesMelIII$Chrom)){
  pMin<-clinesMelIII$pMin[i]
  pMax<-clinesMelIII$pMax[i]
  width<-clinesMelIII$width[i]
  center<-clinesMelIII$center[i]
  deltaM<-clinesMelIII$deltaM[i]
  tauM<-clinesMelIII$tauM[i]
  curve(pMax + (pMin - pMax) * ifelse(deltaM < -(x - center), 1/(1 + 
                                  exp(4 * deltaM/width)) * exp(tauM/(1 + exp(-4 * deltaM/width)) *
                                  ((x - center) * 4/width + 4 * deltaM/width)), ifelse(deltaM < 
                                  x - center, 1 - 1/(1 + exp(4 * deltaM/width)) * exp(-(tauM/(1 +                                                              
				  deltaM/width)), 1/(1 + exp(-((x - center) * 4/width))))),
        	add=T,col=rgb(0.5,0.5,0.5,alpha = 0.5,maxColorValue = 1))
}

# wntA centre
abline(v=49.8  ,col=rgb(204,204,0,max=255))

# Optix centre
abline(v=28.9,col=rgb(192,32,42,max=255))

# Highlight optix and wntA clinesMel:
optix<-clinesMel[clinesMel$Chrom=="Hmel218003o"&clinesMel$Pos==923741,]
curve(optix$pMax + (optix$pMin - optix$pMax) * (1/(1 + exp(-((x - optix$center) * 4/optix$width)))),
      add=T,col=rgb(192,32,42,max=255),lwd=2)

wntA<-clinesMel[clinesMel$Chrom=="Hmel210001o"&clinesMel$Pos==3321549,]
curve(wntA$pMax + (wntA$pMin - wntA$pMax) * (1/(1 + exp(-((x - wntA$center) * 4/wntA$width)))),
      add=T,col=rgb(204,204,0,max=255),lwd=2)

ro<-clinesMel[clinesMel$Chrom=="Hmel213001o"&clinesMel$Pos==10266255,]
curve(ro$pMax + (ro$pMin - ro$pMax) * (1/(1 + exp(-((x - ro$center) * 4/ro$width)))),
      add=T,col=rgb(0,0,204,max=255),lwd=2)

cortex<-clinesMel[clinesMel$Chrom=="Hmel215003o"&clinesMel$Pos==1480717,]
curve(cortex$pMax + (cortex$pMin - cortex$pMax) * (1/(1 + exp(-((x - cortex$center) * 4/cortex$width)))),
      add=T,col=rgb(204,51,255,max=255),lwd=2)
legend("topright",lwd=2,col=c(rgb(192,32,42,max=255),rgb(204,204,0,max=255),
                              rgb(0,0,204,max=255),rgb(204,51,255,max=255)),
                              legend=c(expression(italic("optix")),
                                       expression(italic("WntA")),
                                       expression(italic("Ro")),
                                       expression(italic("cortex"))))


# Plot the distributions of cline widths and centres:
par(mfrow=c(1,2))

hist(c(clinesMel$width,clinesMelIII$width),main="cline widths",breaks=50,xlab="width")
# optix width:
abline(v=15.2,col=rgb(192,32,42,max=255))
# wntA width:
abline(v=24.8 ,col=rgb(204,204,0,max=255))


hist(c(clinesMel$center,clinesMelIII$center),main="cline centres",breaks=50,xlab="centre")
# wntA centre
abline(v=49.8 ,col=rgb(204,204,0,max=255))
# Optix centre
abline(v=28.9,col=rgb(192,32,42,max=255))






