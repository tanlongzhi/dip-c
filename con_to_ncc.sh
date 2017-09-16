#!/bin/bash

input_file="$1"
output_file=${input_file/.con.gz/.ncc}

gunzip -c "$input_file" | awk -F'[\t,]' -v OFS='\t' 'BEGIN {OFS = " "; haplotypes[0] = "pat"; haplotypes[1] = "mat"} {print $1"("haplotypes[$3]")",$2,$2,$2,$2,"+",$4"("haplotypes[$6]")",$5,$5,$5,$5,"+",NR,0,0;}' > "$output_file"