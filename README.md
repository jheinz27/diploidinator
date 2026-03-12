# The Diplinator

Most aligners were not designed for diploid assemblies (e.g. HG002), so when aligning reads to a diploid assembly, the mapping quality for reads may be lower, as there are multiple locations the read can align to well. We have developed a tool to report the best haploid alignment for each read based on the weighted alignment score for all primary and supplemental records (see details below). This tool works for SAM/BAM/CRAM format files or PAF files.

## Installation

Recommended: Download precompiled binary:
```bash
# Todo
```

Build from source:
```bash
git clone https://github.com/jheinz27/diplinator.git
cd diplinator
cargo build --release
./target/release/diplinator
```

## Usage

```
Diplinator: Choose the best alignment to each haploid of a diploid assembly

Usage: diplinator [OPTIONS] <ASM1> <ASM2>

Arguments:
  <ASM1>  asm1 alignment file (sam/bam/cram/paf)
  <ASM2>  asm2 alignment file (sam/bam/cram/paf)

Options:
      --ref1 <FILE>      reference FASTA for cram file (asm1)
      --ref2 <FILE>      reference FASTA for cram file (asm2)
  -1, --s1 <NAME>        label for asm1 sample (used in output file names and summary) [default: asm1]
  -2, --s2 <NAME>        label for asm2 sample (used in output file names and summary) [default: asm2]
  -u, --unmapped <DEST>  where to write reads unmapped in both assemblies: asm1, asm2, or discard [default: asm1] [possible values: asm1, asm2, discard]
      --paf              input files are PAF
      --ms               use ms:i: tag rather than AS:i: for alignment score
  -b, --both             write reads with equal alignment scores to both output files
  -t, --threads <INT>    Total thread pool size (min 4). Multiples of 8 recommended for optimal read/write balance. [default: 8]
  -h, --help             Print help
  -V, --version          Print version
```

## Example Workflow

Diplinator only works on name-sorted files, which is the default [minimap2](https://github.com/lh3/minimap2) output. If you need to use it on coordinate-sorted BAM files, name-sort them first:

```bash
samtools sort -n -o out_name_sort.bam in_index_sort.bam
```

### Diploid assembly alignment

```bash
# If needed, split diploid genome assembly FASTA into respective haplotypes
./separate_haps_fasta hg002v1.1.fa  # writes hg002v1.1.hap1.fa and hg002v1.1.hap2.fa

# Align reads to each haplotype
minimap2 -ax map-hifi -o asm1_alignments.sam hg002v1.1.hap1.fa reads.fastq
minimap2 -ax map-hifi -o asm2_alignments.sam hg002v1.1.hap2.fa reads.fastq

# Run diplinator
diplinator -1 hap1 -2 hap2 asm1_alignments.sam asm2_alignments.sam
# Output: diplinator_hap1.sam  diplinator_hap2.sam

# Merge best alignments into one file (if desired)
samtools merge -@ 12 merged.sam diplinator_hap1.sam diplinator_hap2.sam

# Save as sorted BAM
samtools sort -@ 12 -o merged.bam merged.sam
```

### Comparing different reference genomes

Diplinator can also be used to select best alignments between different reference genomes (e.g. GRCh38 and CHM13):

```bash
diplinator -1 grch38 -2 chm13 grch38_alignments.sam chm13_alignments.sam
# Output: diplinator_grch38.sam  diplinator_chm13.sam
```

### Merging output files

Output files can be merged into one file using samtools:

```bash
samtools merge -@ 12 merged.bam diplinator_asm1.bam diplinator_asm2.bam
```

### CRAM input files

If input files are CRAM format, the original reference genome must be provided:

```bash
diplinator --ref1 asm1_hap.fasta --ref2 asm2_hap.fasta asm1_alignments.cram asm2_alignments.cram
```

## Example PAF Usage

**NOTE:** It is important to use the `--secondary=no --paf-no-hit` flags when aligning with minimap2. If a SAM file is converted to a PAF file with `paftools.js sam2paf`, it will NOT have the required AS tag.

```bash
minimap2 -cx splice -uf -k14 -t 16 --secondary=no --paf-no-hit hg002v1.1.asm1.fasta reads.fastq > out_asm1.paf
minimap2 -cx splice -uf -k14 -t 16 --secondary=no --paf-no-hit hg002v1.1.asm2.fasta reads.fastq > out_asm2.paf

diplinator --paf out_asm1.paf out_asm2.paf
# Output: diplinator_asm1.paf  diplinator_asm2.paf
```

## Scoring Mechanism

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

If the alignment score of the read is equal in both reference genomes, then the "better" alignment is assigned deterministically using a hash of the read name (or written to both output files if `--both` is used).
