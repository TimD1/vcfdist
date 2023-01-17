#!/bin/bash

source includes.sh

# get list of all beds
beds=""
for i in `seq 0 $((${#truth_names[@]}-1))`; do
    beds+=" $data/${truth_names[i]}/$ref_id/${truth_beds[i]}"
done
beds="${beds:1}"


cat $beds | sort -k1,1 -k2,2n -k3,3n -u | bedtools merge -i stdin > all.bed
