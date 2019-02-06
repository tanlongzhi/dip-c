#!/bin/bash

input_file="$1"
output_file=${input_file/.3dg/.dip-c.3dg}


backbone_unit=$(grep "^#unit" "$input_file" | cut -d' ' -f2)

echo "input: distance unit = "${backbone_unit}

grep -v "^#" "$input_file" | sed 's/a	/(pat)	/g; s/b	/(mat)	/g; s/X	/X(mat)	/g; s/Y	/Y(pat)	/g' | awk -F'\t' -v OFS='\t' -v backbone_unit=${backbone_unit} 'BEGIN {norm_factor=1.0/backbone_unit} {print $1,$2,$3*norm_factor,$4*norm_factor,$5*norm_factor}' > "$output_file"

