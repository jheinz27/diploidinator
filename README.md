
# The diplinator
Most aligners were not designed for diploid assemblies (eg. HG002), so when aligning reads to a diploid assembly, the mapping quality for reads may be lower, as there are multiple locations the read can align to well. We have developed a tool to report the best haploid aligment for each read based on the weighted alignment score for all primary and supplemental records (see details below). This tool works for SAM/BAM/CRAM format files or PAF files.


## Installation
Recommended: Download precompiled binary:
```
Todo
```

Build from source:
```
git clone https://github.com/jheinz27/diploidinator.git
cd diploidinator
cargo build --release
./target/release/diploidinator
```

## Usage
```
Diploidinator: Choose the best alignment to each haploid of a diploid assembly for every read

Usage: diploidinator [OPTIONS] <ASM1> <ASM2>

Arguments:
  <ASM1>  asm1 alignment file (sam/bam/cram/paf)
  <ASM2>  asm2 alignment file (sam/bam/cram/paf)

Options:
      --ref1 <FILE>      reference FASTA for cram file (asm1)
      --ref2 <FILE>      reference FASTA for cram file (asm2)
      --s1 <NAME>        label for asm1 sample (used in output file names and summary) [default: asm1]
      --s2 <NAME>        label for asm2 sample (used in output file names and summary) [default: asm2]
  -o, --out <PREFIX>     prefix of output files [default: diploidinator_out]
      --paf              input files are PAF
  -t, --threads <INT>    Total thread pool size (min 4). Multiples of 8 recommended for optimal read/write balance. [default: 8]
  -h, --help             Print help
  -V, --version          Print version
```



## Example Workflow
The Diploidinator will only work on a name sorted file, which is the default [minimap2](https://github.com/lh3/minimap2) output. However, if it is desired to use the diploidinator on index sorted bam files, it will be necessary to name sort them with `samtools sort -n -o out_name_sort.bam in_index_sort.bam`

```
#if needed, split diploid genome assembly fasta into respective haplotypes
./separate_haps_fasta hg002v1.1.fa #writes hg002v1.1.hap1.fa and hg002v1.1.hap2.fa
#align reads to each haplotype
minimap2 -ax map-hifi -o asm1_alignments.sam hg002v1.1.hap1.fa reads.fastq
minimap2 -ax map-hifi -o asm2_alignments.sam hg002v1.1.hap2.fa reads.fastq
#run diplinator 
./diploidinator -o diplinator_out asm1_alignments.sam asm2_alignments.sam

merge best aligments into one file (if desired)

//save as sorted bam 
samtools sort
samtools merge -@ 12 merged.sam diplinator_hap1.sam diplinator_hap2.sam
```

With custom sample labels (e.g. grch38 and chm13):
```
./diploidinator --s1 grch38 --s2 chm13 -o diplinator_out grch38_alignments.sam chm13_alignments.sam
# Output: diplinator_out_grch38.bam  diplinator_out_chm13.bam
```

### Merging output files

diploidinator output files can be merged into one file with:
```
samtools merge -@ 12 merged.bam diplinator_out_asm1.bam diplinator_out_asm2.bam
```

### CRAM input files
if input files are CRAM format, the original referece genome will need to be provided:
```
./diploidinator --ref1 asm1_hap.fasta --ref2 asm2_hap.fasta -o diplinator_out asm1_alignments.cram asm2_alignments.cram
```



## Example PAF Usage
NOTE: It is important to use the `--secondary=no --paf-no-hit` flags when aligning with minimap2. If a SAM file is converted to a PAF file with `paftools.js sam2paf` it will NOT have the required AS tag.
```
minimap2 -cx splice -uf -k14 -t 16 --secondary=no --paf-no-hit hg002v1.1.asm1.fasta read.fastq > out_asm1.paf
minimap2 -cx splice -uf -k14 -t 16 --secondary=no --paf-no-hit hg002v1.1.asm2.fasta reads.fastq > out_asm2.paf
./diploidinator --paf out_asm1.paf out_asm2.paf
```

## Scoring mechanism used to comparing alignments

Let $n$ be the number of primary and supplementary alignments for a given read to that reference genome.

Let $I_i = [r_i^{\mathrm{start}}, r_i^{\mathrm{end}})$ denote the interval on the read covered by alignment $i$.

Let

$$
B = \left| \bigcup_{i=0}^{n} I_i \right|
$$

be the number of **read bases** covered by at least one alignment.

Let $a_i$ be the alignment score for alignment $i$.

Let $l_i$ be the alignment length in read coordinates for alignment $i$.

$$
S = B \cdot \frac{\sum_{i=0}^{n} a_i}{\sum_{i=0}^{n} l_i}
$$

If the alignment score of the read is equal in both reference genomes, then the "better" alignment is assigned randomly by taking the last digit of the hash of the read name.
