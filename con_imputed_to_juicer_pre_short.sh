#!/bin/bash

input_file="$1"
output_file=${input_file/.con.gz/.juicer.txt.gz}

gunzip -c "$input_file" | awk -F'[\t,]' -v OFS='\t' '{hom1=$1 "(" ($3==0?"pat":"mat") ")"; hom2=$4 "(" ($6==0?"pat":"mat") ")"; if(hom1<hom2){print 0,hom1,$2,0,0,hom2,$5,1;}else{print 0,hom2,$5,0,0,hom1,$2,1;}}' | sort -k2,2 -k6,6 | gzip -c > "$output_file"