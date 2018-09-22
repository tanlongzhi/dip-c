#!/bin/bash

input_file="$1"
output_file=${input_file/.3dg/.dip-c.3dg}

grep -v "^#" "$input_file" | sed 's/a	/(pat)	/g; s/b	/(mat)	/g; s/X	/X(mat)	/g; s/Y	/Y(pat)	/g' > "$output_file"
