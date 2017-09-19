#!/bin/bash

sample_folder="1CDX1-413"
snp_file="snps/cast_129.txt.gz"
metac_path="."
preset_name="mf"

${metac_path}/metac seg -v ${snp_file} ${sample_folder}/srt.bam | gzip -c > ${sample_folder}/phased.seg.gz
${metac_path}/metac con ${sample_folder}/phased.seg.gz | gzip -c > ${sample_folder}/raw.con.gz
${metac_path}/metac dedup ${sample_folder}/raw.con.gz | gzip -c > ${sample_folder}/dedup.con.gz
${metac_path}/metac reg -p ${preset_name} ${sample_folder}/dedup.con.gz | gzip -c > ${sample_folder}/seg.con.gz
