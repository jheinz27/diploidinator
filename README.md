# The Diplinator

Most aligners were not designed for diploid assemblies (e.g. HG002), so when aligning reads to a diploid assembly, the mapping quality for reads may be lower, as there are multiple locations the read can align to well. We have developed a tool to report the best haploid alignment for each read based on a weighted alignment score computed across all primary and supplemental records (see details below). For each read, Diplinator also reports a haplotype assignment quality (HAPQ) score that captures how confidently the read can be assigned to one haplotype over the other. This tool works for SAM/BAM/CRAM format files or PAF files.

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
  -1, --s1 <NAME>          label for asm1 sample (used in output file names and summary) [default: asm1]
  -2, --s2 <NAME>          label for asm2 sample (used in output file names and summary) [default: asm2]
      --paf                input files are PAF
      --ms                 use ms:i: tag rather than AS:i: for alignment score
  -b, --both               write reads with equal alignment scores to both output files
  -u, --unmapped <DEST>    where to write reads unmapped in both assemblies: asm1, asm2, or discard [default: asm1] [possible values: asm1, asm2, discard]
      --ref1 <FILE>        reference FASTA for cram file (asm1)
      --ref2 <FILE>        reference FASTA for cram file (asm2)
      --match-sc <FLOAT>   per-base match score from aligner scoring scheme (e.g. minimap2 default is 2.0 for long reads) [default: 2.0]
      --no-hapq            skip HAPQ score calculation and hq tag output (e.g. for comparing grch38 vs chm13)
  -t, --threads <INT>      Total thread pool size (min 4). Multiples of 8 recommended for optimal read/write balance. [default: 8]
  -h, --help               Print help
  -V, --version            Print version
```

Each output record is annotated with an `hq:i:` tag carrying the HAPQ score (see [HAPQ](#hapq-haplotype-assignment-quality)), unless `--no-hapq` is set.

## Example Workflow

Diplinator only works on name-sorted files, which is the default [minimap2](https://github.com/lh3/minimap2) output. If you need to use it on coordinate-sorted BAM files, name-sort them first:

```bash
samtools sort -n -o out_name_sort.bam in_index_sort.bam
```

### Diploid assembly alignment

```bash
# If needed, split diploid genome assembly FASTA into respective haplotypes
separate_haps_fasta -1 MATERNAL -2 PATERNAL hg002v1.1.fa
# writes hg002v1.1.MATERNAL.fa and hg002v1.1.PATERNAL.fa

# Align reads to each haplotype
minimap2 -ax map-hifi -o asm1_alignments.sam hg002v1.1.MATERNAL.fa reads.fastq
minimap2 -ax map-hifi -o asm2_alignments.sam hg002v1.1.PATERNAL.fa reads.fastq

# Run diplinator
diplinator -1 mat -2 pat asm1_alignments.sam asm2_alignments.sam
# Output: diplinator_mat.sam  diplinator_pat.sam

# Merge best alignments into one file (if desired)
samtools merge -@ 12 merged.sam diplinator_mat.sam diplinator_pat.sam

# Save as sorted BAM
samtools sort -@ 12 -o merged.bam merged.sam
```

### Comparing different reference genomes

Diplinator can also be used to select best alignments between different reference genomes (e.g. GRCh38 and CHM13). For this use case, the HAPQ score is generally not meaningful, so pass `--no-hapq` to skip its calculation:

```bash
diplinator --no-hapq -1 grch38 -2 chm13 grch38_alignments.sam chm13_alignments.sam
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

For each read, Diplinator computes a single weighted alignment score per assembly using all primary and supplementary alignments (secondary alignments are passed through to the output but ignored when scoring).

Let $n$ be the number of primary and supplementary alignments for a given read to that reference genome.

Let $L$ be the full read length (sum of query-consuming CIGAR operations on a non-secondary record, or the `qlen` field of the PAF record).

Let $I_i = [r_i^{\mathrm{start}}, r_i^{\mathrm{end}})$ denote the interval on the read covered by alignment $i$.

Let

$$
B = \left| \bigcup_{i=1}^{n} I_i \right|
$$

be the number of **read bases** covered by at least one alignment.

Let $a_i$ be the alignment score for alignment $i$ (`AS:i:` by default, or `ms:i:` if `--ms` is set).

Let $l_i$ be the alignment length in read coordinates for alignment $i$.

$$
S = \frac{\sum_{i=1}^{n} a_i}{\sum_{i=1}^{n} l_i} \cdot B \cdot \frac{B}{L}
$$

The first factor is the average alignment score per aligned base (across all non-secondary segments). It is multiplied by the number of unique read bases covered $B$ and then scaled by the read coverage fraction $B/L$, so that reads which align over a large fraction of their length are weighted more heavily than reads which align only over a small portion.

For each read, the assembly with the higher $S$ wins; its full alignment cluster (including secondary alignments) is written to the corresponding output file. If $S$ is equal in both assemblies, the "better" assignment is determined by a hash of the read name (deterministic and reproducible), or the read is written to both output files when `--both` is used.

## HAPQ (haplotype assignment quality)

For each read assigned to a winning haplotype, Diplinator reports a HAPQ score in the `hq:i:` tag of the output record. HAPQ is a Phred-like confidence (0-60) that the read was assigned to the correct haplotype. The calculation is modeled on BWA-MEM's `mem_approx_mapq_se`.

Let $S_w$ and $S_l$ be the weighted alignment scores of the winning and losing assemblies, $m$ the per-base match score (`--match-sc`, default `2.0`), and $k$ the number of non-secondary alignments (splits) on the winning side.

$$
\Delta = \frac{S_w - S_l}{m}
\qquad
p_\text{split} = \begin{cases} 1 & k \le 3 \\ \dfrac{3}{k} & k > 3 \end{cases}
$$

$$
\text{HAPQ} = \mathrm{clamp}\bigl(6.02 \cdot \Delta \cdot p_\text{split},\ 0,\ 60\bigr)
$$

$\Delta$ is approximately the difference, in matching bases, between the winning and losing alignments. The split penalty $p_\text{split}$ down-weights reads with more than three split alignments, which often fall in complex/repetitive regions where haplotype assignment is less reliable.

Special cases:
- Read mapped in only one assembly: HAPQ = 60.
- Read tied between assemblies (winner = `Both`): HAPQ = 0.
- Read unmapped in both assemblies: no `hq` tag is written.

If `--no-hapq` is set, HAPQ is not computed and no `hq:i:` tag is added (recommended when the two inputs are not haplotypes of the same sample, e.g. GRCh38 vs CHM13).
