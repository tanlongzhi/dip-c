# Test Data Generation

## Source

All reads were extracted from **GSM3271358** (Tan et al., *Science* 2018), a single GM12878 cell (Cell 11) sequenced with Dip-C. The raw reads (SRR7226708, 30.5M paired-end read pairs) were mapped to **hs37d5** (GRCh37) using:

```bash
bwa mem -5SP -t 64 \
  -R "@RG\tID:SRR7226708\tPL:ILLUMINA\tSM:SRR7226708" \
  /path/to/hs37d5/bwa_index/genome.fa \
  SRR7226708_1.fastq SRR7226708_2.fastq \
| samtools view -bS - \
| samtools sort -o SRR7226708.bam -
```

This matches the standard Dip-C mapping pipeline: bwa 0.7.17 with `-5SP` flags (Hi-C mode), producing a coordinate-sorted BAM. Note that `-5SP` means **0% of reads are marked as properly paired** — all paired reads have `is_proper_pair=False`.

## Extraction process

A Python script (`extract_test_reads.py`) using pysam scanned chr22 of the full BAM to find reads matching specific criteria, then extracted **all alignments** for those read names from the entire BAM (including mates and supplementary alignments on other chromosomes). The script performed a full linear scan of all 66.5M alignments to ensure no alignments were missed.

## test_seg.bam / test_seg.bam.bai

Coordinate-sorted, indexed BAM containing **91 alignments** from **26 unique read names** (26 read1 + 26 read2 = 52 primary + 39 supplementary).

### Read categories

**Core chimeric reads (SA tag present, segments on chr22):**

| Read name | Category | Description |
|---|---|---|
| SRR7226708.27522991 | chimeric_inter_chr | SA points to chr14 (inter-chromosomal contact) |
| SRR7226708.6662870 | chimeric_inter_chr | SA points to chr1 |
| SRR7226708.28041940 | chimeric_intra_far | SA on chr22, >10Mb apart (intra-chromosomal, far) |
| SRR7226708.20299677 | chimeric_intra_far | SA on chr22, >10Mb apart |
| SRR7226708.9938350 | chimeric_multi_sa | 2+ SA entries (3+ total segments per mate) |
| SRR7226708.13122872 | chimeric_multi_sa | 2+ SA entries |
| SRR7226708.10101251 | chimeric_reverse_sa | SA alignment on reverse strand |
| SRR7226708.15009667 | chimeric_reverse_sa | SA alignment on reverse strand |

**Discordant read pairs (no SA tag, paired, not proper pair):**

| Read name | Category | Description |
|---|---|---|
| SRR7226708.5948288 | discordant | Both mates mapped, not proper pair, no SA tag |
| SRR7226708.7251007 | discordant | " |
| SRR7226708.10232484 | discordant | " |
| SRR7226708.10681329 | discordant | " |

Note: With `-5SP` mapping, these are NOT "negative controls" — seg processes them as discordant reads (since `is_proper_pair=False` for all reads). However, they are cleaned out during `seg_data.clean()` because they produce only 1 segment per mate and both mates enter via pass 1, but the resulting reads still have < 2 qualifying segments after quality filtering.

**Quality filter edge cases:**

| Read name | Category | Description |
|---|---|---|
| SRR7226708.2363794 | low_mapq_chimeric | Has SA tag but primary MAPQ < 20; seg filters out its segments |
| SRR7226708.11247069 | low_mapq_chimeric | Same; also appears as a supplementary alignment |
| SRR7226708.30346024 | high_nm | NM/(aligned length) > 0.05; seg filters via `max_nm_per_bp` |
| SRR7226708.20582799 | high_nm | Same |
| SRR7226708.10043407 | chimeric_sa_low_mapq | Primary passes filters (MAPQ=60), but one SA entry has low MAPQ — only the bad SA segment is dropped |
| SRR7226708.25556901 | chimeric_sa_low_mapq | Same |

**Both mates chimeric:**

| Read name | Category | Description |
|---|---|---|
| SRR7226708.9079525 | both_mates_chimeric | Both R1 and R2 have SA tags; tests no double-counting |
| SRR7226708.12799820 | both_mates_chimeric | Same |

**Supplementary alignments (flag 0x800):**

| Read name | Category | Description |
|---|---|---|
| SRR7226708.19035906 | supplementary | Supplementary alignment on chr22; seg skips these in pass 1 but catches the primary via its SA tag |
| SRR7226708.11247069 | supplementary | Overlaps with low_mapq_chimeric (same read) |

**Duplicate-flagged read (flag 0x400):**

One read was synthetically marked as duplicate during extraction (flag 0x400 set on all its alignments). This read had no SA tag and was not otherwise selected. Seg skips duplicate-flagged reads entirely. Visible in flagstat as "2 duplicates" (both mates of the same read).

**SNP phasing reads:**

| Read name | Category | Description |
|---|---|---|
| SRR7226708.8572487 | snp_paternal | Overlaps SNP at 22:16059753; base matches paternal allele (T), bq=40 |
| SRR7226708.3067002 | snp_paternal | Overlaps SNP at 22:16065774; base matches paternal allele (A), bq=39 |
| SRR7226708.988747 | snp_maternal + snp_low_baseq | Overlaps SNP at 22:16059753; base matches maternal allele (A), but bq=2 (should be filtered by min_baseq=20) |
| SRR7226708.1164589 | snp_maternal | Overlaps SNP at 22:16059753; base matches maternal allele (A), bq=39 |
| SRR7226708.2707682 | snp_neither | Overlaps SNP at 22:16066157; base is G but alleles are T/C (neither match, should be skipped) |

Note: These SNP reads don't appear in the seg output because they lack SA tags and their discordant pair segments don't survive `seg_data.clean()`. They are included in the BAM for potential future tests of the pileup/phasing logic in isolation.

### Flagstat summary

```
91 + 0 in total
52 + 0 primary
39 + 0 supplementary
2 + 0 duplicates
90 + 0 mapped (98.90%)
51 + 0 primary mapped (98.08%)
0 + 0 properly paired (0.00%)
1 + 0 singletons (1.92%)
10 + 0 with mate mapped to a different chr
```

## test_snps.txt

Tab-delimited SNP file (9 entries), format: `chr  pos(1-based)  paternal_base  maternal_base`.

These are **not** from the NA12878 SNP file — only 1 real NA12878 SNP on chr22 overlapped the 12 reads that survive seg filtering (the reads are 250bp and there are ~31K SNPs across ~51Mb of chr22, so overlap is sparse). Instead, these SNPs were constructed by examining the actual bases in the reads at specific positions using pileup, then setting the paternal/maternal alleles to create controlled test scenarios:

| Position | Pat | Mat | Read(s) hit | Base in read | bq | Expected behavior |
|---|---|---|---|---|---|---|
| 22:16208451 | A | G | .6662870 (bq=27), .10101251 (bq=16) | A | 27, 16 | .6662870 phased paternal; .10101251 filtered (bq < 20) |
| 22:16900471 | G | A | .9938350, .13122872, .20299677 | G | 40 | All phased paternal |
| 22:17282541 | A | G | .10043407 R1 | A | 39 | Phased paternal |
| 22:17282851 | C | A | .10043407 R2 | A | 39 | Phased maternal |
| 22:17519051 | T | G | .25556901 R2 | G | 39 | Phased maternal |
| 22:28740821 | A | G | .9938350, .13122872, .20299677 | G | 39-40 | All phased maternal |
| 22:32393651 | T | C | .9079525 R1 | C | 40 | Phased maternal |
| 22:32617834 | C | T | .9079525 R1+R2 | C | 37-38 | Phased paternal (this is also a real NA12878 SNP) |
| 22:37896551 | C | G | .12799820 | A | 38-39 | Neither allele matches — skipped |

## test_seg.seg

Output of `dip-c seg test_seg.bam` (no phasing). Contains **12 reads** with **44 segments**, 0% phased.

The 12 reads that survive are all chimeric (have SA tags with qualifying segments). The 4 reads removed by `seg_data.clean()` had < 2 segments after quality filtering (low MAPQ or high NM reads where all segments were filtered out).

## test_seg_phased.seg

Output of `dip-c seg -v test_snps.txt test_seg.bam`. Same 12 reads, 44 segments, **31.82% phased** (14/44 segments).

Phased segments are marked with `0` (paternal) or `1` (maternal) in the last field of each segment. Unphased segments remain `.`.

## test_seg.con

Output of `dip-c con test_seg_phased.seg`. Contains **18 contacts**, 44.44% legs phased. Contacts include both inter-chromosomal (chr1-22, chr2-22, chr14-22, chr18-22, chr22-hs37d5) and intra-chromosomal (chr22-22) entries with phased haplotype information.

## Python 3 compatibility notes

The original dip-c code was written for Python 2. Running with Python 3 requires these fixes:

1. **Dict mutation during iteration** (`classes.py`): `SegData.clean()` and multiple `ConData`/`G3dData` methods delete dict entries while iterating `.keys()`. In Python 3, `.keys()` returns a view, not a copy. Fix: wrap in `list()` or pre-collect keys to delete.

2. **Binary vs text file mode** (`seg.py`, `con.py`): SNP and seg files opened with `"rb"` produce bytes in Python 3, causing silent failures in string comparisons (e.g., `b"T" != "T"`). Fix: use `"rt"` / `"r"` mode.

3. **Pileup stepper** (`seg.py`): In pysam >=0.22, the default pileup stepper (`samtools`) sets `ignore_orphans=True`, which filters reads whose mates aren't nearby. Dip-C reads are often "orphans" (chimeric mates on different chromosomes). Fix: use `stepper='all', min_base_quality=0` (seg does its own base quality filtering).
