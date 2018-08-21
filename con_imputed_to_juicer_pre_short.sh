#!/bin/bash

input_file="$1"
output_file=${input_file/.con.gz/.juicer.txt.gz}

gunzip -c "$input_file" | awk -F'[\t,]' -v OFS='\t' '{print 0,$1 "(" ($3==0?"pat":"mat") ")",$2,0,0,$4 "(" ($6==0?"pat":"mat") ")",$5,1}' | sort -k2,2 -k6,6 | gzip -c > "$output_file"