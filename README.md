
# The diploid-inator 


## Installation 

## Usage
```
Choose the best alignment to each haploid of a diploid assembly

Usage: diploidinator_sams [OPTIONS] --mat <FILE> --pat <FILE>

Options:
  -m, --mat <FILE>      hap1.sam/bam/cram
  -p, --pat <FILE>      hap2.sam/bam/cram
      --ref-mat <FILE>  reference FASTA for cram file
      --ref-pat <FILE>  reference FASTA for cram file
  -o, --out <PREFIX>    hap2.sam/bam/cram [default: diploidinator_out]
      --paf             input files are PAF
  -t, --threads <INT>   Number of threads to use for BAM file decompression [default: 4]
  -h, --help            Print help
  -V, --version         Print version
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

## example usage 
