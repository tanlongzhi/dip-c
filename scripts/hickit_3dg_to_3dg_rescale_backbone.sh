#!/bin/bash

input_file="$1"
output_file=${input_file/.3dg/.dip-c.3dg}


backbone_unit=$(grep -v "^#" "$input_file" | awk '{chr=$1;loc=$2;x=$3;y=$4;z=$5;if(NR>=2){if(chr==prev_chr){print sqrt((prev_x-x)^2+(prev_y-y)^2+(prev_z-z)^2)}};prev_chr=chr;prev_loc=loc;prev_x=x;prev_y=y;prev_z=z;}' | sort -n | awk '{dist[NR]=$0}END{if(NR%2){print dist[(NR+1)/2]}else{print (dist[(NR/2)]+dist[(NR/2)+1])/2.0}}')

echo "input: median backbone distance = "${backbone_unit}

grep -v "^#" "$input_file" | sed 's/a	/(pat)	/g; s/b	/(mat)	/g; s/X	/X(mat)	/g; s/Y	/Y(pat)	/g' | awk -F'\t' -v OFS='\t' -v backbone_unit=${backbone_unit} 'BEGIN {norm_factor=2.0/backbone_unit} {print $1,$2,$3*norm_factor,$4*norm_factor,$5*norm_factor}' > "$output_file"

