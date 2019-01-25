#!/bin/bash

input_file="$1"
output_file=${input_file/.cen/.cen.leg}

awk -v OFS=',' '{print $1,0,0; print $1,$3,0; print $1,$2,0; print $1,0,1; print $1,$3,1; print $1,$2,1;}' ${input_file} > ${output_file}