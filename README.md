# diploidinator

out_mat=${b/.fastq.gz/_to_hg002.MAT.paf}
out_pat=${b/.fastq.gz/_to_hg002.PAT.paf}
out_haps=${b/.fastq.gz/_to_hg002_hap_merge.paf}
echo $out_mat
minimap2 -cx splice -uf -k14 -t 16 --secondary=no --paf-no-hit /hlilab/jakob/genomes/hg002v1.1.MATERNAL.fasta $fq > $out_mat
echo $out_pat
minimap2 -cx splice -uf -k14 -t 16 --secondary=no --paf-no-hit /hlilab/jakob/genomes/hg002v1.1.PATERNAL.fasta $fq > $out_pat
python /hlilab/jakob/cell_lines_pairs/manuscript_v1/diploidinator_v2.py $out_mat $out_pat > $out_haps
