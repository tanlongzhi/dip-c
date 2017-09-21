#!/bin/bash

sample_folder="1CDX1-413"
snp_file="snps/cast_129.txt.gz"
metac_path="."
preset_name="mf"

${metac_path}/metac seg -v ${snp_file} ${sample_folder}/srt.bam | gzip -c > ${sample_folder}/phased.seg.gz
${metac_path}/metac con ${sample_folder}/phased.seg.gz | gzip -c > ${sample_folder}/raw.con.gz
${metac_path}/metac dedup ${sample_folder}/raw.con.gz | gzip -c > ${sample_folder}/dedup.con.gz
${metac_path}/metac reg -p ${preset_name} ${sample_folder}/dedup.con.gz | gzip -c > ${sample_folder}/reg.con.gz
${metac_path}/metac clean -c2 ${sample_folder}/reg.con.gz | gzip -c > ${sample_folder}/clean.con.gz
${metac_path}/metac impute ${sample_folder}/clean.con.gz | gzip -c > ${sample_folder}/impute.con.gz

${metac_path}/con_to_ncc.sh ${sample_folder}/impute.con.gz
${metac_path}/nuc_dynamics.sh ${sample_folder}/impute.ncc 0.1
${metac_path}/metac impute3 -3 ${sample_folder}/impute.3dg ${sample_folder}/clean.con.gz | gzip -c > ${sample_folder}/impute3.round1.con.gz

${metac_path}/con_to_ncc.sh ${sample_folder}/impute3.round1.con.gz
${metac_path}/nuc_dynamics.sh ${sample_folder}/impute3.round1.ncc 0.1
${metac_path}/metac impute3 -3 ${sample_folder}/impute3.round1.3dg ${sample_folder}/clean.con.gz | gzip -c > ${sample_folder}/impute3.round2.con.gz

${metac_path}/con_to_ncc.sh ${sample_folder}/impute3.round2.con.gz
${metac_path}/nuc_dynamics.sh ${sample_folder}/impute3.round2.ncc 0.1
${metac_path}/metac impute3 -3 ${sample_folder}/impute3.round2.3dg ${sample_folder}/clean.con.gz | gzip -c > ${sample_folder}/impute3.round3.con.gz

${metac_path}/metac clean3 -c ${sample_folder}/impute3.round2.con.gz ${sample_folder}/impute3.round2.3dg > ${sample_folder}/impute3.round2.clean.3dg