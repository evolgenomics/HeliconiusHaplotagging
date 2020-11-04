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
nice realSFS fst stats2 melpomene.fst.idx -win 10000 -step 2500 -type 0 > melpomene.folded.fst.10kb.2.5kb.idx
 
# Liftover the melpomene window positions to erato reference genome positions
file=melpomene.folded.fst.10kb.2.5kb.idx
liftOver <(awk 'NR>1{print $2,$3,$3+1000,$2"_"$3"_"$5}' $file) \
  hmel2.5-helera1_demo.largeGap.14k.filtered.chain.gz \
  $file.eratoLiftOver.tmp $file.eratoLiftOver.unmapped -multiple

# Add a header and remove additional matches (index>1)
sed -e '1iscaffold\tpos\tpos2\tHmelScaff\tMelPos\tFst\tindex' \
 -e 's/_/\t/g' $file.eratoLiftOver.tmp | awk '{if($7<2 || $7=="index") print}'  > $file.eratoLiftOver
rm $file.eratoLiftOver.tmp # remove the temporary file

