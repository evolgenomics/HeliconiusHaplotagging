#!/bin/bash
file=$1
grep BLOCK -A 1 $file| awk '/BLOCK/ {span = $9;frags=$11;phased=$7};/Herato/ {print $4"\t"$5"\t"$5+span"\t"$4":"$5"-"$5+span":"phased":"span":"frags}' > ${file/.output/.bed}
