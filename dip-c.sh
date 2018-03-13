#!/bin/bash

# paths
sample_folder="gm12878_03"
snp_file="snps/gm12878.txt.gz"
dip-c_path="."

# parameters
preset_name="hf"
clean_count="5"
impute_clean_count="2"
impute3_clean_count="2"

# preparation
${dip-c_path}/dip-c seg -v ${snp_file} ${sample_folder}/aln.bam | gzip -c > ${sample_folder}/phased.seg.gz
${dip-c_path}/dip-c con ${sample_folder}/phased.seg.gz | gzip -c > ${sample_folder}/raw.con.gz
${dip-c_path}/dip-c dedup ${sample_folder}/raw.con.gz | gzip -c > ${sample_folder}/dedup.con.gz
${dip-c_path}/dip-c reg -p ${preset_name} ${sample_folder}/dedup.con.gz | gzip -c > ${sample_folder}/reg.con.gz
${dip-c_path}/dip-c clean -c${clean_count} ${sample_folder}/reg.con.gz | gzip -c > ${sample_folder}/clean.con.gz
${dip-c_path}/dip-c impute -C${impute_clean_count} ${sample_folder}/clean.con.gz | gzip -c > ${sample_folder}/impute.con.gz

# 3d reconstruction
${dip-c_path}/con_to_ncc.sh ${sample_folder}/impute.con.gz
${dip-c_path}/nuc_dynamics.sh ${sample_folder}/impute.ncc 0.1
${dip-c_path}/dip-c impute3 -C${impute3_clean_count} -p ${preset_name} -3 ${sample_folder}/impute.3dg ${sample_folder}/clean.con.gz | gzip -c > ${sample_folder}/impute3.round1.con.gz
${dip-c_path}/dip-c clean3 -c ${sample_folder}/impute.con.gz ${sample_folder}/impute.3dg > ${sample_folder}/impute.clean.3dg

${dip-c_path}/con_to_ncc.sh ${sample_folder}/impute3.round1.con.gz
${dip-c_path}/nuc_dynamics.sh ${sample_folder}/impute3.round1.ncc 0.1
${dip-c_path}/dip-c impute3 -C${impute3_clean_count} -p ${preset_name} -3 ${sample_folder}/impute3.round1.3dg ${sample_folder}/clean.con.gz | gzip -c > ${sample_folder}/impute3.round2.con.gz
${dip-c_path}/dip-c clean3 -c ${sample_folder}/impute3.round1.con.gz ${sample_folder}/impute3.round1.3dg > ${sample_folder}/impute3.round1.clean.3dg

${dip-c_path}/con_to_ncc.sh ${sample_folder}/impute3.round2.con.gz
${dip-c_path}/nuc_dynamics.sh ${sample_folder}/impute3.round2.ncc 0.1
${dip-c_path}/dip-c impute3 -C${impute3_clean_count} -p ${preset_name} -3 ${sample_folder}/impute3.round2.3dg ${sample_folder}/clean.con.gz | gzip -c > ${sample_folder}/impute3.round3.con.gz
${dip-c_path}/dip-c clean3 -c ${sample_folder}/impute3.round2.con.gz ${sample_folder}/impute3.round2.3dg > ${sample_folder}/impute3.round2.clean.3dg

${dip-c_path}/con_to_ncc.sh ${sample_folder}/impute3.round3.con.gz
${dip-c_path}/nuc_dynamics.sh ${sample_folder}/impute3.round3.ncc 0.02
${dip-c_path}/dip-c impute3 -C${impute3_clean_count} -p ${preset_name} -3 ${sample_folder}/impute3.round3.3dg ${sample_folder}/clean.con.gz | gzip -c > ${sample_folder}/impute3.round4.con.gz
${dip-c_path}/dip-c clean3 -c ${sample_folder}/impute3.round3.con.gz ${sample_folder}/impute3.round3.3dg > ${sample_folder}/impute3.round3.clean.3dg

${dip-c_path}/con_to_ncc.sh ${sample_folder}/impute3.round4.con.gz
${dip-c_path}/nuc_dynamics.sh ${sample_folder}/impute3.round4.ncc 0.02
${dip-c_path}/dip-c clean3 -c ${sample_folder}/impute3.round4.con.gz ${sample_folder}/impute3.round4.3dg > ${sample_folder}/impute3.round4.clean.3dg