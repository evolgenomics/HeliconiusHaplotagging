# Bash code to compute 10kb window averages of FST with angsd using the STITCH/HAPCUT2 genotype likelihoods for Heliconius melpomene and liftover to erato positions

# Run each scaffold:
for SCAFF in `sed -n -e "$SLURM_ARRAY_TASK_ID p" scaffolds`  # note, this was run as an array script on a cluster. The $SLURM_ARRAY_TASK_ID is just the job index (i.e. 1, 2, 3...)
do
 FILE=run164_merged_${SCAFF}.PL.AD.HAPCUT2.inf0.05.corr

 # Extract the malleti and the plesseni samples from the vcf file
 bcftools view -Oz -o ${FILE}.malleti.vcf.gz -S malleti.samples ../vcf/$FILE.vcf.gz
 bcftools view -Oz -o ${FILE}.plesseni.vcf.gz -S plesseni.samples ../vcf/$FILE.vcf.gz

 # Compute the 1D-SFS fore ach subspecies
 angsd -vcf-gl ${FILE}.malleti.vcf.gz -anc $ref -out $FILE.malleti -dosaf 1 
 angsd -vcf-gl ${FILE}.plesseni.vcf.gz -anc $ref -out $FILE.plesseni -dosaf 1

done

# Concatenate the 1D-SFS of the different scaffolds
realSFS cat -outnames malleti run164_merged_*malleti.saf.idx  
realSFS cat -outnames plesseni run164_merged_*plesseni.saf.idx

# Make FOLDED 2d-SFS
realSFS plesseni.merged.saf.idx malleti.merged.saf.idx -P 8 -fold 1 > plesseni-malleti.folded.ml

# Compute overall FST
nice realSFS fst index plesseni.merged.saf.idx malleti.merged.saf.idx \
 -fstout melpomene -whichFST 1 -sfs plesseni-malleti.folded.ml  

# Compute sliding window Fst in 10kb windows
nice realSFS fst stats2 melpomene.fst.idx -P 6 -win 10000 -step 10000 -type 0 > melpomene.folded.fst.10kb.idx

# Run a Hidden Markov Model (HMM) to define regions of elevated differentiation

# Compute z scores in R
R
fst<-read.table("fst.10kb.idx",header=T)
zscore <- function(x) {
  x<- as.numeric(as.character(x))
  (x-mean(x))/sd(x) }
z<-zscore(fst$Fst)
write.table(z,file="zscores.10kb.txt",row.names=F,col.names=F,quote=F)

# Add a header
sed -i '1ix' zscores.10kb.txt

# Run a HMM over these z scores specifying 3 cores
# Transition and emission probabilities are optimized with Baum-Welch algorithm
# ML sequence of state is inferred with the Viterbi algorithm
# Note, I fixed the means of the state distributions to 0 and 3 and the sd to 1
Rscript HMM_zscores_2norm_fixedMeanSd.R zscores.10kb.txt 3

# Add the HMM states to the 10 kb Fst values
paste fst.10kb.idx <(cut -d" "  -f 2 zscores.10kb_2state_HMMstates.txt) > fst.10kb.idx.withHMM

# liftover the window positions to H. erato reference genome positions
file=fst.10kb.idx.withHMM
cp $file ../eratoPos
cd ../eratoPos
liftOver <(awk 'NR>1{print $2,$3,$3+1000,$2"_"$3"_"$5"_"$6}' $file) \
  hmel2.5-helera1_demo.largeGap.14k.filtered.chain.gz \
  $file.eratoLiftOver.tmp $file.eratoLiftOver.unmapped -multiple
 
sed -e '1iscaffold\tpos\tpos2\tHmelScaff\tMelPos\tFst\tHMM\tindex' \
 -e 's/_/\t/g' $file.eratoLiftOver.tmp | awk '{if($8<2 || $8=="index") print}'  > $file.eratoLiftOver
rm $file.eratoLiftOver.tmp

