#!/bin/bash

input_file="$1"
output_file=${input_file/.pairs.gz/.con.gz}

gunzip -c "$input_file" | grep -v "^#" | awk -v OFS='\t' '{chr1=substr($2,1,length($2)-1);hap1=substr($2,length($2),1);if(hap1=="X"){chr1=$2;hap1=1}else if(hap1=="Y"){chr1=$2;hap1=0}else{hap1=hap1=="b"};chr2=substr($4,1,length($4)-1);hap2=substr($4,length($4),1);if(hap2=="X"){chr2=$4;hap2=1}else if(hap2=="Y"){chr2=$4;hap2=0}else{hap2=hap2=="b"};print chr1","$3","hap1, chr2","$5","hap2;}' | gzip -c > "$output_file"