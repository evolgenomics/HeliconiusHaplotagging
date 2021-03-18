#!/bin/bash
tag=$1
nHaps=$3
nSnps=$2
rho=$4
reps=$5

for rep in `seq $reps`; do mspms --recombination $rho $nHaps $nSnps $rep --precision 7 -T > simARG/mspms_${tag}_rep${rep}.out; done
