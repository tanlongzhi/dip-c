#!/bin/bash

input_file="$1"
output_file=${input_file/.len/.hom.len}

awk -v OFS='\t' '{print $1"(pat)",$2;print $1"(mat)",$2;}' ${input_file} > ${output_file}