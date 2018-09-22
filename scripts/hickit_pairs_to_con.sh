#!/bin/bash

input_file="$1"
output_file=${input_file/.pairs.gz/.con.gz}

gunzip -c "$input_file" | grep -v "^#" | awk -v OFS='\t' '{print $2","$3","$8, $4","$5","$9}' | gzip -c > "$output_file"
