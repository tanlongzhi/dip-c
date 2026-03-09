# Workflow

## Typical Workflow

In our latest work, both the main Dip-C algorithm and 3D modeling are now carried out with [hickit](https://github.com/lh3/hickit), a much faster and more careful implementation. Below is a typical workflow of such combined use of [hickit](https://github.com/lh3/hickit) and this repo (with mm10 as an example genome):

```sh
# align reads
bwa mem -5SP genome.fa R1.fq.gz R2.fq.gz | gzip > aln.sam.gz # for Nextera
#seqtk mergepe R1.fq.gz R2.fq.gz | pre-meta - | bwa mem -5SP -p genome.fa - | gzip > aln.sam.gz # for META

# extract segments
hickit.js sam2seg -v snp.txt.gz aln.sam.gz | hickit.js chronly -y - | gzip > contacts.seg.gz # for female
#hickit.js sam2seg -v snp.txt.gz aln.sam.gz | hickit.js chronly - | hickit.js bedflt par.bed - | gzip > contacts.seg.gz # for male

# resolve haplotypes via imputation
hickit -i contacts.seg.gz -o - | bgzip > contacts.pairs.gz
hickit -i contacts.pairs.gz -u -o - | bgzip > impute.pairs.gz

# generate 3D structures (with 3 replicates)
for rep in `seq 1 3`
do
  hickit -s${rep} -M -i impute.pairs.gz -Sr1m -c1 -r10m -c2 -b4m -b1m -O 1m.${rep}.3dg -b200k -O 200k.${rep}.3dg -D5 -b50k -O 50k.${rep}.3dg -D5 -b20k -O 20k.${rep}.3dg
done

# convert from hickit to dip-c formats, and remove repetitive regions from 3D structures
scripts/hickit_pairs_to_con.sh contacts.pairs.gz
scripts/hickit_impute_pairs_to_con.sh impute.pairs.gz
for rep in `seq 1 3`
do
  scripts/hickit_3dg_to_3dg_rescale_unit.sh 20k.${rep}.3dg
  dip-c clean3 -c impute.con.gz 20k.${rep}.dip-c.3dg > 20k.${rep}.clean.3dg # remove repetitive (contact-less) regions
done

# align replicate structures and calculate RMSD (overall value in .log file)
dip-c align -o aligned.20k. 20k.[1-3].clean.3dg 2> 20k.align.log > 20k.align.color

# convert to juicebox format for interactive viewing
# raw contacts
java -Xmx2g -jar juicer_tools.jar pre -n contacts.pairs.gz contacts.hic mm10
# haplotype-resolved contacts
scripts/con_imputed_to_juicer_pre_short.sh impute.con.gz
java -Xmx2g -jar juicer_tools.jar pre -n impute.juicer.txt.gz impute.hic color/mm10.chr.hom.len

# calculate single-cell chromatin compartment values along the genome
dip-c color2 -b1000000 -H -c color/mm10.cpg.1m.txt -s contacts.con.gz > cpg_b1m.color2 # contact-based
dip-c color -c color/mm10.cpg.20k.txt -s3 20k.1.clean.3dg > cpg_s3.color # 3D-structure-based

# calculate radial positioning
dip-c color -C 20k.1.clean.3dg > C.color

# color by chromosome number and visualize as mmCIF (viewable with pymol)
dip-c color -n color/mm10.chr.txt 20k.1.clean.3dg | dip-c vis -c /dev/stdin 20k.1.clean.3dg > 20k.1.clean.n.cif
```

## Stand-alone Workflow (Old)

The original workflow, which was before the development of [hickit](https://github.com/lh3/hickit) and uses a slightly modified version of [nuc_dynamics](https://github.com/TheLaueLab/nuc_dynamics) for 3D modeling, can be found in an archived [document](https://github.com/tanlongzhi/dip-c/blob/master/README_old.md).
