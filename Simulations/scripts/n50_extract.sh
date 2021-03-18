#!/bin/bash
#usage - ./n50_extract.sh <input hapcut2_to_bed.bed> [n threshold]
file=$1
if [ -n "$2" ]; then num=$2; else num=50; fi
sum=`sed 's/:/\t/g' $file | sort -k 7,7nr | /fml/chones/local/bin/datamash sum 7`
threshold=`echo $sum | awk '{print sprintf("%.2f", $1*'$num'/100)}'`
sed 's/:/\t/g' $file | sort -k 7,7nr | awk 'NR ==1  && $7 >= '$threshold'{print $1"\t"$2"\t"$3"\t"$4":"$5"\t"$6"\t"$7"\t"$8} ; {sum+=$7};sum<'$threshold'{print $1"\t"$2"\t"$3"\t"$4":"$5"\t"$6"\t"$7"\t"$8}' | sort -k 1,1 -k 2,2n | tee ${file/.bed/.n$num.bed}; 
