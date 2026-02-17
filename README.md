
# The diploid-inator 
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

Usage: diploidinator_sams [OPTIONS] --mat <FILE> --pat <FILE>

Options:
  -m, --mat <FILE>      hap1.sam/bam/cram
  -p, --pat <FILE>      hap2.sam/bam/cram
      --ref-mat <FILE>  reference FASTA for cram file (mat)
      --ref-pat <FILE>  reference FASTA for cram file (pat)
  -o, --out <PREFIX>    prefix of output files [default: diploidinator_out]
      --paf             input files are PAF
  -t, --threads <INT>   Total thread pool size (min 4). Multiples of 8 recommended for optimal read/write balance. [default: 8]
  -h, --help            Print help
  -V, --version         Print version
```



## Example Usage 
The Diploidinator will only work on a name sorted file, which is the default [minimap2](https://github.com/lh3/minimap2) output. However, if it is desired to use the diploidinator on index sorted bam files, it will be necessary to name sort them with `samtools sort -n -o out_name_sort.bam in_index_sort.bam`

```
minimap2 -ax map-hifi -o mat_alignments.sam genomes/hg002v1.1.MATERNAL.fasta reads.fastq
minimap2 -ax map-hifi -o pat_alignments.sam genomes/hg002v1.1.PATERNAL.fasta reads.fastq
./diploidinator -m mat_alignments.sam -p pat_alignments.sam -o diplinator_out
```
### Merging output files

diploidinator output files can be merged into one file with: 
```
samtools merge -@ 12 merged.bam diplinator_out_mat.bam diplinator_out_pat.bam
```

### CRAM input files 
if input files are CRAM format, the original referece genome will need to be provided: 
```
./diploidinator -m mat_alignments.sam -p pat_alignments.sam -o diplinator_out --ref-mat maternal_hap.fasta --ref-pat paternal_hap.fasta 
```



## Example PAF Usage
NOTE: It is important to use the `--secondary=no --paf-no-hit` flags when aligning with minimap2. If a SAM file is converted to a PAF file with `paftools.js sam2paf` it will NOT have the required AS tag. 
```
minimap2 -cx splice -uf -k14 -t 16 --secondary=no --paf-no-hit hg002v1.1.MATERNAL.fasta read.fastq > out_mat.paf
minimap2 -cx splice -uf -k14 -t 16 --secondary=no --paf-no-hit hg002v1.1.PATERNAL.fasta reads.fastq > out_pat.paf 
./diploidinator -m out_mat.paf -m out_pat.paf --paf
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
