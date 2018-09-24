#!/bin/bash

input_file="$1"
output_file=${input_file/.3dg/.dip-c.3dg}

num_beads=$(grep -v "^#" "$input_file" | wc -l | tr -d ' ')

grep -v "^#" "$input_file" | sed 's/a	/(pat)	/g; s/b	/(mat)	/g; s/X	/X(mat)	/g; s/Y	/Y(pat)	/g' | awk -F'\t' -v OFS='\t' -v num_beads=${num_beads} 'BEGIN {norm_factor=num_beads^(1/3)/10} {print $1,$2,$3*norm_factor,$4*norm_factor,$5*norm_factor}' > "$output_file"

