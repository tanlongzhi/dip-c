# File Formats

## Phased SNPs

To resolve the two haplotypes, a list of phased single-nucleotide polymorphisms (SNPs) is required in a tab-delimited format (chromosome, coordinate, paternal nucleotide, maternal nucleotide).

An example SNP file for the GM12878 cell line (in hg19 coordinates), based on the 2016-1.0 release of [Illumina Platinum Genomes](https://www.illumina.com/platinumgenomes.html), is provided as `snps/NA12878.txt.gz`. The first few lines are shown below:

```
1	534324	G	T
1	550515	T	C
1	565419	G	C
1	655642	G	A
1	666543	G	A
1	676118	T	C
1	705882	A	G
1	727529	A	G
1	732772	G	A
1	742825	A	G
```

## Read Segments

The basic output of a chromatin conformation capture (3C/Hi-C) experiment is sequencing reads (or read pairs) containing more than one genomic segments. In each of these reads (sometimes known as *chimeric* reads), genomic segments far away in the linear genome (or even from different chromosomes) are joined by proximity-based ligation.

These segments are recorded in a `.seg` file as an intermediate format. Each line represents a read (or read pair) in a tab-delimited format (read name, segment 1, segment 2 ...).

Each segment is recorded as a comma-delimited string: `.` for read 1 and `m` (short for "mate") for read 2, start coordinate in the read, end coordinate in the read, chromosome, start coordinate in the genome, end coordinate in the genome, strand of the genome (`+` or `-`), haplotype (`.` for unknown, `0` for paternal, and `1` for maternal). A segment has a known haplotype if it carries one or more phased SNPs.

An example `.seg` file:

```
HWI-D00433:595:HLYW7BCXY:1:1209:15116:100489    .,0,211,5,94770308,94770519,+,. m,57,211,5,97167374,97167528,-,0        m,0,61,5,94770541,94770602,-,.
HWI-D00433:595:HLYW7BCXY:1:1113:20520:48210     .,0,193,21,23683758,23683951,-,0        .,188,255,21,25124149,25124216,-,.
HWI-D00433:595:HLYW7BCXY:2:2207:6749:75115      .,39,92,15,71782216,71782269,-,.        .,0,43,15,62676681,62676724,+,.
HWI-D00433:595:HLYW7BCXY:1:1103:9994:47614      .,44,278,1,19116987,19117221,+,.        .,0,48,1,19108099,19108147,-,.
HWI-D00433:595:HLYW7BCXY:2:1216:14071:49584     .,0,114,11,120680807,120680921,+,.      .,109,211,11,120689618,120689720,-,.    m,0,211,11,120689228,120689439,+,.
```

## Contact Legs

A "leg" is an endpoint of a read segment — it forms one of the two legs of a chromatin contact. More generally, a leg can be any single point in the genome.

Each leg is recorded as a comma-delimited string: chromosome, coordinate, haplotype.

An example `.leg` file:

```
1,948359,.
1,1192624,.
1,2561820,.
1,2836242,.
1,2954969,1
1,3114198,.
1,3343831,.
1,3455767,.
1,3518062,.
1,3540154,1
```

## Contacts

Chromatin contacts are a crucial concept in 3C/Hi-C. A contact refers to the proximity-based ligation of two genomic segments far away in the linear genome (or even from different chromosomes). A contact is defined as the two adjoining endpoints (legs) of two different segments in a same read (or read pair).

Each contact is recorded as a tab-delimited line: leg 1, leg 2.

An example `.con` file:

```
1,858641,.      1,861338,.
1,861471,.      1,862872,.
1,918037,.      1,1024147,.
1,918249,1      1,1231502,.
1,921617,0      1,956928,.
1,922873,.      1,926783,.
1,923319,.      1,957711,.
1,946196,.      1,1235547,.
1,948480,.      1,1133615,.
1,959161,.      1,962343,.
```

!!! note "Subtleties"
    - When a read contains more than two segments, contacts are defined between **all** pairs of segments (similar to bulk studies of multi-way contacts).
    - The coordinate of a leg is exact if it resides in the middle of read 1 or read 2. It is approximate if it resides in the unread gap between read 1 and read 2.
    - The two legs of a contact are interchangeable. To avoid ambiguity, leg 1 < leg 2.
    - Directionalities are ignored. [hickit](https://github.com/lh3/hickit) preserves this additional information using the [`.pairs` format](https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md).

## 3D Genomes

The primary output of the Dip-C algorithm is the 3D structure of a single-cell genome. Following the definition in [nuc_dynamics](https://github.com/TheLaueLab/nuc_dynamics), the 3D structure of each chromosome is represented by 3D coordinates of regularly spaced points (0 kb, 10 kb, 20 kb, 30 kb ...) along the chromosome. 3D coordinates of points elsewhere will be linearly interpolated from the given points.

Each 3D genome is recorded in a tab-delimited format (chromosome with `(pat)` for paternal and `(mat)` for maternal, coordinate, x, y, z).

An example `.3dg` file:

```
1(mat)  1420000 0.791377837067  10.9947291355   -13.1882897693
1(mat)  1440000 -0.268241283699 10.5200875887   -13.0896257278
1(mat)  1460000 -1.3853075236   10.5513787498   -13.1440142173
1(mat)  1480000 -1.55984101733  11.4340829129   -13.6026301209
1(mat)  1500000 -0.770991778399 11.4758488546   -14.5881137222
1(mat)  1520000 -0.0848245107875        12.2624690808   -14.354289628
1(mat)  1540000 -0.458643807046 12.5985791771   -13.4701149287
1(mat)  1560000 -0.810322906201 12.2461643989   -12.3172933413
1(mat)  1580000 -2.08211172035  12.8886838656   -12.8742007778
1(mat)  1600000 -3.52093948201  13.1850935438   -12.4118684428
```

## Genomic Regions

A `.reg` file performs a similar role to a [BED file](https://genome.ucsc.edu/FAQ/FAQformat.html#format1), but with haplotype information. This format can be used to exclude regions of copy-number (CN) gains or losses of heterozygosity (LOHs) from a `.con` file, to set haplotypes in regions of CN losses in a `.con` file, or to extract regions of interest from a `.3dg` file.

Each region is recorded as a tab-delimited line: chromosome, haplotype, start coordinate (`.` for the start of the chromosome), end coordinate (`.` for the end of the chromosome).

An example `.reg` file:

```
1	.	.	.
2	1	232800000	.
2	.	238500000	.
4	.	.	.
5	.	165800000	167800000
5	1	167800000	.
6	1	.	.
11	.	.	.
16	.	40000000	.
19	.	.	.
```
