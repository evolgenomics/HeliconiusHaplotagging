#!/bin/bash
file=$1
awk '/chr/ && $2 == 0 {print $4"\t"$5-1"\t"$5"\t"$6"\t"$7"\t0/0\t1/1"};/chr/ && $2 == 1 {print $4"\t"$5-1"\t"$5"\t"$7"\t"$6"\t0/0\t1/1"}' $file > ${file/.output/.phaseSet}
