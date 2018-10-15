#!/bin/bash

input_file="$1"
output_file=${input_file/.con.gz/.juicer.txt.gz}

gunzip -c "$input_file" | awk -F'[\t,]' -v OFS='\t' '{print 0,$1,$2,0,0,$4,$5,1}' | gzip -c > "$output_file"