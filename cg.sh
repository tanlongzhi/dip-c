#!/bin/bash
#SBATCH -n 1
#SBATCH -t 3-00:00
#SBATCH -p general
#SBATCH --mem-per-cpu=8000

genome_file="$1"
window_size="$2"

half_window_size=$(expr ${window_size} / 2)

tmpdir=`mktemp -d --tmpdir=.`

bedtools makewindows -g ${genome_file}.fai -w ${window_size} -s ${half_window_size} | awk '($3%'${window_size}'=='${half_window_size}'){print $0;}' > ${tmpdir}/windows.bed
while IFS='' read -r line || [[ -n "$line" ]]; do
    printf '%s' "${line}" | bedtools getfasta -fi ${genome_file} -bed stdin | grep -v "^>" | fold -1 | awk '{if($0=="A"||$0=="T"){atsum+=1;}else if($0=="C"||$0=="G"){cgsum+=1;};} END {if(atsum+cgsum>0){print cgsum/(atsum+cgsum);}else{print -1;}}'
done < ${tmpdir}/windows.bed > ${tmpdir}/values.txt

paste ${tmpdir}/windows.bed ${tmpdir}/values.txt | awk -v OFS='\t' '($4>=0){print $1,($2+$3)/2,$4}'

rm -r ${tmpdir}
