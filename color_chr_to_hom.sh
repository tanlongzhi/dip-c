#!/bin/bash

input_file="$1"
output_file=${input_file/.txt/.hom.txt}

awk -v OFS='\t' '{print $1"(pat)",$2,$3;print $1"(mat)",$2,$3;}' ${input_file} > ${output_file}