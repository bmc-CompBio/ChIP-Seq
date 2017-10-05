#!/bin/sh
# workflow.sh

# Iterate through files
for FILENAME in *.sra
do
  # extract file basename
	FILEBASE=`echo ${FILENAME} | sed -e "s/.sra//g"`
  
  # convert sra to fastq
	fastq-dump ${FILEBASE}.sra
  
  # run quality control
	fastqc ${FILEBASE}.fastq

  # run bowtie2 alignment
	bowtie2  -x Drosophila_melanogaster.BDGP6.dna.selected -U ${FILEBASE}.fastq -p 4 -S ${FILEBASE}.sam 2> ${FILEBASE}.stats 
  
  # convert sam to bam, sort by coordinate and index
	samtools view -bS -q 0 ${FILEBASE}.sam | samtools sort - | tee ${FILEBASE}.bam | samtools index - ${FILEBASE}.bam.bai
  
  # convert bam to bedgraph
	genomeCoverageBed -ibam ${FILEBASE}.bam  -bg -g Drosophila_melanogaster.BDGP6.dna.selected.chromInfo.txt > ${FILEBASE}.bedgraph
  
  # convert bedgraph to tdf
	igvtools toTDF -z 5 ${FILEBASE}.bedgraph ${FILEBASE}.tdf Drosophila_melanogaster.BDGP6.dna.selected.fa

done