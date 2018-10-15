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
    printf '%s' "${line}" | awk -F $'\t' 'BEGIN {OFS = FS} {print $1,$2,$3+1;}' | bedtools getfasta -fi ${genome_file} -bed stdin | grep -v "^>" | awk '{for(i=1;i<length($0);i++){subseq=substr($0,i,2); split(subseq,subseqchars,"");if((subseqchars[1]!="A"&&subseqchars[1]!="C"&&subseqchars[1]!="G"&&subseqchars[1]!="T")||(subseqchars[2]!="A"&&subseqchars[2]!="C"&&subseqchars[2]!="G"&&subseqchars[2]!="T")){continue;};sum+=1;if(subseq=="CG"){cgsum+=1;}}} END {if(sum>0){print cgsum/sum;}else{print -1;}}'
done < ${tmpdir}/windows.bed > ${tmpdir}/values.txt

paste ${tmpdir}/windows.bed ${tmpdir}/values.txt | awk -v OFS='\t' '($4>=0){print $1,($2+$3)/2,$4}'

rm -r ${tmpdir}
