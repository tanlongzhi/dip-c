#!/bin/bash

resolution_array=(8 4 2 0.4 0.2 0.1 0.04 0.02 0.01)

input_ncc_file="$1"
final_resolution="$2"
nuc_dynamics_path="/Users/tanlongzhi/Desktop/Zones/HiC/Code/nuc_dynamics-dip-harmonicbackbone"

# make array string
resolution_string=""
for resolution in "${resolution_array[@]}"
do
    resolution_string=$(echo ${resolution_string} ${resolution})
    if [ "${resolution}" = "${final_resolution}" ]
    then
        break
    fi
done
echo "simulation steps (Mb):" ${resolution_string}

${nuc_dynamics_path}/nuc_dynamics ${input_ncc_file} -o ${input_ncc_file/.ncc/.3dg} -temps 20 -s ${resolution_string}