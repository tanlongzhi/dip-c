#!/bin/bash

input_file="$1"
output_file=${input_file/.con.gz/.con.txt}

gunzip -c "$input_file" | awk -F'[\t,]' -v OFS='\t' '{print $1,$2,$3,$4,$5,$6}' | sed 's/X/20/g; s/Y/21/g; s/\./-1/g; s/chr//g' > "$output_file"