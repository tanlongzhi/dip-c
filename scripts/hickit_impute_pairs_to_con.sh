#!/bin/bash

input_file="$1"
output_file=${input_file/.pairs.gz/.con.gz}

threshold="0.75"

gunzip -c "$input_file" | grep -v "^#" | awk -v OFS='\t' '{if($10>='${threshold}'){hap1=0;hap2=0;}else if($11>='${threshold}'){hap1=0;hap2=1;}else if($12>='${threshold}'){hap1=1;hap2=0;}else if($13>='${threshold}'){hap1=1;hap2=1;}else{next;}; print $2","$3","hap1, $4","$5","hap2;}' | gzip -c > "$output_file"