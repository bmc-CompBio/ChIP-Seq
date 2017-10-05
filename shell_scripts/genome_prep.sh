#!/bin/sh
# genome_prep.sh 

# Download genome
wget ftp://ftp.ensembl.org/pub/release-90/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.dna.toplevel.fa.gz

# Extract fasta file
gunzip Drosophila_melanogaster.BDGP6.dna.toplevel.fa.gz

# Select chromosomes
samtools faidx Drosophila_melanogaster.BDGP6.dna.toplevel.fa X 2L 2R 3L 3R > Drosophila_melanogaster.BDGP6.dna.selected.fa

# Index new fasta file
samtools faidx Drosophila_melanogaster.BDGP6.dna.selected.fa

# generate chromInfo table with chromosome length
cut -f1,2 Drosophila_melanogaster.BDGP6.dna.selected.fa.fai > Drosophila_melanogaster.BDGP6.dna.selected.chromInfo.txt

# build bowtie2 index
bowtie2-build Drosophila_melanogaster.BDGP6.dna.selected.fa Drosophila_melanogaster.BDGP6.dna.selected
