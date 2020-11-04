#!Rscript --vanilla

####################################################################################
# R-script that runs an HMM on z-scores with two normally distributed hidden states
# The mean and standard deviations are fixed
# sd=empirical sd for both states, means=empirical mean and 90th quantile
# high zscores represent high FST values
# This script was written by Joana Meier on June-2017 and is based on a script by David Marques
#
####################################################################

#########
# USAGE #
#########

# Rscript HMM_zscores_2norm_fixedMeanSd.R <infile> <ncores>
# infile contains z-scores with the header "x"
# Read input parameters
args<-commandArgs(trailingOnly=T)

# Check if 2 arguments were given
if(length(args)!=2){
  stop("Usage: Rscript HMM_zscores_2norm_fixedMeanSd.R <filewithzscores.txt> <ncores>\n\t\t- filewithzscores.txt is file with probabilities between 0 and 1 (quantile/q-values from outlier analysis), that has a header 'x'\n\t\t- ncores is the number of cores allowed for computation") 
}

##################
# Load libraries #
##################

library("HiddenMarkov")
library("foreach")
library("doParallel")

################################
# Load input data (z-scores)   #
################################
data<-read.table(file=args[1], header=T, quote="\"")
zscore<-data$x

print(paste("Read",length(zscore),"sites"))

# Define base name for output files
elements<-strsplit(args[1],split="\\.")
base<-paste(elements[[1]][-length(elements[[1]])],collapse=".")
rm(elements)


#######################################################################################
# Define Mstep for HMM (syntax follows the rules by the Hidden Markov R package)
# This function modifies the BaumWelch parameter estimation by fixing the means and sd
# Note: if dthmm is run with "nsd", it expects the functions Mstep.nsd,
#       pnsd (distribution function), dnsd (density) and rnsd (random number generator)
# Here pnsd, dnsd and rnsd are set to default
#######################################################################################

# Define custom Maximisation step, derived from normal distribution
# with fixed standard deviation and means of both states
# leave the mean and standard deviations as initially set
Mstep.nsd <- function(x, cond, pm, pn){
  #   this function is a modified version of Mstep.norm
  #   here the mean is fixed to the values specified in pm$mean
  nms <- sort(names(pm))
  n <- length(x)
  m <- ncol(cond$u)
  if (all(nms==c("mean", "sd"))){
    mean <- pm$mean
    sd <- pm$sd
    return(list(mean=mean, sd=sd))
  }
}
# Set the distribution function, density and random number generator to default (normal distribution)
rnsd <- rnorm
dnsd <- dnorm
pnsd <- pnorm
qnsd <- qnorm


####################################
# parallelize parameter estimation #
# - parameter estimation from 1000 #
#   random starting parameters     #
####################################

# Define the number of cores to be used
registerDoParallel(cores=args[2])

# Run parameter estimation in parallel on args[2] cores
print("Performing parameter estimation with 1000 iterations")
parout<-foreach(i=1:1000,.packages="HiddenMarkov",.combine='c') %dopar% {
  print(i);set.seed(i)
  # Sample random initial parameters for the transition matrix (trans),
  #   marginal/initial probabilities (init) and for the state means
  prob<-runif(3,min=0,max=1)
  prob2<-runif(1,min=0,max=(1-prob[3]))
  trans<-matrix(c(prob[1],1-prob[1],
                  prob[2],(1-prob[2])),byrow=T,nrow=2)
  init<-c(prob[3],prob2,(1-(prob[3]+prob2)))

  # Define the fixed standard deviations and means:
  #sdevs<-c(sd(zscore),sd(zscore))
  #means<-c(mean(zscore),quantile(zscore,probs = 0.1))
 
  # With theoretical sd and means:
  sdevs<-c(1,1) 
  means<-c(0,3) # fixing the means of the state distributions

  # Build Hidden Markov Model with random intial and fixed parameters
  myhmm<-dthmm(zscore,trans,init,"nsd",list(mean=means,sd=sdevs),discrete=F)

  # Run Baum-Welch algorithm, 3 runs to find maximal estimates
  # Baum-Welch configuration
  a<-bwcontrol(maxiter=1000,tol=1e-07,prt=F,posdiff=F)
  bwhmm2<-try(BaumWelch(myhmm,control=a),silent=T)
  bwhmm1<-try(BaumWelch(bwhmm2,control=a),silent=T)
  bwhmm<-try(BaumWelch(bwhmm1,control=a),silent=T)
  
  # Output parameters
  if(length(bwhmm)>1){
    c(bwhmm$Pi,bwhmm$delta,bwhmm$pm$mean,bwhmm$pm$sd,bwhmm$LL,bwhmm$iter)
  }else{
    rep(NA,12)
  }
}

hmmpar<-as.data.frame(matrix(parout,ncol=12,byrow=T))
names(hmmpar)<-c("a11","a21","a12","a22","i1","i2","m1","m2","sd1","sd2","LL","iter")
write.table(hmmpar,file=paste(base,"_2state_HMMparest",sep=""))
hmmpar<-read.table(paste(base,"_2state_HMMparest",sep=""),stringsAsFactors=F)

########################################
## EVALUATION OF PARAMETER ESTIMATION ##
########################################

# Plot parameter estimates
#pdf(file=paste(base,"_2state_parest_plot.pdf",sep=""),height=16,width=10)
#par(mfrow=c(7,3))
#attach(hmmpar)
#for(k in colnames(hmmpar)){
#  plot(hmmpar$LL,get(k),xlab="Likelihood",ylab=k,pch=20,col="#88888888",xlim=c(min(hmmpar$LL,na.rm=T),max(hmmpar$LL,na.rm=T)))
#  abline(h=with(hmmpar,get(k)[LL==max(LL,na.rm=T)]),col="red")
#}
#dev.off()

# Store the best initial parameters (the one with the highest likelihood)
bestpar<-head(hmmpar[with(hmmpar,order(-hmmpar$LL)),],1)

print(paste("best paramters: ",bestpar))
print("Running another parameter optimization round to ensure the maximum likelihood is reached")

# Run another parameter optimization round to ensure the maximum likelihood is reached
{
  # Define custom Maximisation step, derived from normal distribution
  # with fixed standard deviations and means
  Mstep.nsd <- function(x, cond, pm, pn){
    #   this function is a modified version of Mstep.norm
    #   here the mean is fixed to the values specified in pm$mean
    nms <- sort(names(pm))
    n <- length(x)
    m <- ncol(cond$u)
    if (all(nms==c("mean", "sd"))){
      mean <- pm$mean
      sd <- pm$sd
      return(list(mean=mean, sd=sd))
    }
  }
  rnsd <- rnorm
  dnsd <- dnorm
  pnsd <- pnorm
  qnsd <- qnorm

  # Build Hidden Markov Model with best parameters
  myhmm<-dthmm(zscore,matrix(c(bestpar$a11,bestpar$a12,
                            bestpar$a21,bestpar$a22),byrow=T,nrow=2),
                c(bestpar$i1,bestpar$i2),"nsd",
                list(mean=c(bestpar$m1,bestpar$m2),sd=c(bestpar$sd1,bestpar$sd2)),discrete=F)
  
  # Run Baum-Welch algorithm 3 times in a row to maximize parameter estimates
  # Baum-Welch configuration
  a<-bwcontrol(maxiter=1000,tol=1e-07,prt=F,posdiff=F)
  bwhmm2<-try(BaumWelch(myhmm,control=a),silent=T)
  bwhmm1<-try(BaumWelch(bwhmm2,control=a),silent=T)
  bwhmm<-try(BaumWelch(bwhmm1,control=a),silent=T)  
}

# Save new best parameter
bestpar<-c(bwhmm$Pi,bwhmm$delta,bwhmm$pm$mean,bwhmm$pm$sd,bwhmm$LL,bwhmm$iter)
names(bestpar)<-c("a11","a21","a12","a22","i1","i2","m1","m2","sd1","sd2","LL","iter")
write.table(bestpar,file=paste(base,"_2state_bestparameters.txt",sep=""))

# Plot the data (z-scores) against the inferred state distributions
pdf(file=paste(base,"_zscore_vs_2statedistributions.pdf",sep=""))
 par(mfrow=c(1,1))
 hist(zscore,breaks=50,col="grey",main="Data vs. Emission & Transition Probabilities",border=F,xlab="z-transformed outlier probability",freq = F)
 for(i in 1:length(bwhmm$pm$mean)){
   points(seq(par("usr")[1],par("usr")[2],length.out=1000),dnorm(seq(par("usr")[1],par("usr")[2],length.out=1000),mean=bwhmm$pm$mean[i],sd=bwhmm$pm$sd[i]),type="l",col=c(4,1)[i],lwd=2)
 }
 legend("topright",legend=c(paste("state1"," mean1=",round(bwhmm$pm$mean[1],2)," sd1=",round(bwhmm$pm$sd[1],2)),
                           paste("state2"," mean2=",round(bwhmm$pm$mean[2],2)," sd2=",round(bwhmm$pm$sd[2],2)),
                           paste("1>1: ",round(bwhmm$Pi[1,1],2),"; 1>2: ",round(bwhmm$Pi[1,2],2)),
                           paste("2>1: ",round(bwhmm$Pi[2,1],2),"; 2>2: ",round(bwhmm$Pi[2,2],2))),
       fill=c(4,1,NA,NA),border=c(1,1,0,0),bty="n")
dev.off()

##############################
## RECONSTRUCTION OF STATES ##
##############################

# Reconstruction of states with the Viterbi algorithm
print("Reconstructing HMM states")
states<-Viterbi(bwhmm)

# Output HMM results
write.table(states,file=paste(base,"_2state_HMMstates.txt",sep=""))
