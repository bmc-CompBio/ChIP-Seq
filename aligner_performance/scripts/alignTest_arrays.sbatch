#! /bin/bash
#
# alignTest.sbatch
#
#SBATCH -J alignTest_array           # A single job name for the array
#SBATCH -p slim18                    # Partition
#SBATCH -n 18                        # 18 cores
#SBATCH -N 1                         # nodes required
#SBATCH -t 0-24:00                   # running time
#SBATCH --mem 50000                  # Memory request
#SBATCH -o alignTest_%A_%a.out       # Standard output
#SBATCH -e alignTest_%A_%a.err       # Standard error

FILENAME=`ls *.fastq | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`

module load ngs/samtools/1.3.1
module load ngs/bowtie1/1.1.2
module load ngs/bowtie2/2.2.9
module load ngs/bwa/0.7.16
module load ngs/STAR/2.5.3a

no_cores=16

## get the file name without extension
FILEID=`echo ${FILENAME} | cut -d\. -f1`

mkdir ${FILEID}_out
cd ${FILEID}_out

### bowtie1_default
method="bowtie1_default"
start_time=`date +%s`
bowtie /work/data/genomes/fly/Drosophila_melanogaster/UCSC/dm6/Sequence/BowtieIndex/genome -S -p ${no_cores} ../${FILEID}.fastq | samtools view -bS -o ${FILEID}.bam
samtools sort -m 2G -@ ${no_cores} -o ${FILEID}.s.bam ${FILEID}.bam
echo run time ${method} is $(expr `date +%s` - $start_time) s
mv ${FILEID}.s.bam ${FILEID}.${method}.bam
samtools index ${FILEID}.${method}.bam
rm ${FILEID}.bam

### bowtie1_unique
method="bowtie1_unique"
start_time=`date +%s`
bowtie /work/data/genomes/fly/Drosophila_melanogaster/UCSC/dm6/Sequence/BowtieIndex/genome -S -m 1 -p ${no_cores} ../${FILEID}.fastq | samtools view -bS -o ${FILEID}.bam
samtools sort -m 2G -@ ${no_cores} -o ${FILEID}.s.bam ${FILEID}.bam
echo run time ${method} is $(expr `date +%s` - $start_time) s
mv ${FILEID}.s.bam ${FILEID}.${method}.bam
samtools index ${FILEID}.${method}.bam
rm ${FILEID}.bam

### bowtie2_default
method="bowtie2_default"
start_time=`date +%s`
bowtie2 -p ${no_cores} -x /work/data/genomes/fly/Drosophila_melanogaster/UCSC/dm6/Sequence/Bowtie2Index/genome -U ../${FILEID}.fastq | samtools view -bS -o ${FILEID}.bam
samtools sort -m 2G -@ ${no_cores} -o ${FILEID}.s.bam ${FILEID}.bam
echo run time ${method} is $(expr `date +%s` - $start_time) s
mv ${FILEID}.s.bam ${FILEID}.${method}.bam
samtools index ${FILEID}.${method}.bam
rm ${FILEID}.bam

### bowtie2_unique
method="bowtie2_unique"
start_time=`date +%s`
bowtie2 -p ${no_cores} -x /work/data/genomes/fly/Drosophila_melanogaster/UCSC/dm6/Sequence/Bowtie2Index/genome -U ../${FILEID}.fastq > ${FILEID}.sam
samtools view -h ${FILEID}.sam | grep -v "XS:i:" | samtools view -b -o ${FILEID}.bam
samtools sort -m 2G -@ ${no_cores} -o ${FILEID}.s.bam ${FILEID}.bam
echo run time ${method} is $(expr `date +%s` - $start_time) s
mv ${FILEID}.s.bam ${FILEID}.${method}.bam
samtools index ${FILEID}.${method}.bam
rm ${FILEID}.bam
rm ${FILEID}.sam

### bowtie2_q0
method="bowtie2_q0"
start_time=`date +%s`
bowtie2 -p ${no_cores} -x /work/data/genomes/fly/Drosophila_melanogaster/UCSC/dm6/Sequence/Bowtie2Index/genome -U ../${FILEID}.fastq | samtools view -bS -q 0 -o ${FILEID}.bam
samtools sort -m 2G -@ ${no_cores} -o ${FILEID}.s.bam ${FILEID}.bam
echo run time ${method} is $(expr `date +%s` - $start_time) s
mv ${FILEID}.s.bam ${FILEID}.${method}.bam
samtools index ${FILEID}.${method}.bam
rm ${FILEID}.bam
rm ${FILEID}.sam

### bwa mem
method="bwa_mem"
start_time=`date +%s`
bwa mem -t ${no_cores} /work/data/genomes/fly/Drosophila_melanogaster/UCSC/dm6/Sequence/BWAIndex/genome.fa ../${FILEID}.fastq | samtools view -bS -o ${FILEID}.bam
samtools sort -m 2G -@ ${no_cores} -o ${FILEID}.s.bam ${FILEID}.bam
echo run time ${method} is $(expr `date +%s` - $start_time) s
mv ${FILEID}.s.bam ${FILEID}.${method}.bam
samtools index ${FILEID}.${method}.bam
rm ${FILEID}.bam

### bwa mem
method="bwa_mem_unique"
start_time=`date +%s`
bwa mem -t ${no_cores} -c 1 /work/data/genomes/fly/Drosophila_melanogaster/UCSC/dm6/Sequence/BWAIndex/genome.fa ../${FILEID}.fastq | samtools view -bS -o ${FILEID}.bam
samtools sort -m 2G -@ ${no_cores} -o ${FILEID}.s.bam ${FILEID}.bam
echo run time ${method} is $(expr `date +%s` - $start_time) s
mv ${FILEID}.s.bam ${FILEID}.${method}.bam
samtools index ${FILEID}.${method}.bam
rm ${FILEID}.bam


### bwa aln
method="bwa_aln"
start_time=`date +%s`
bwa aln -t ${no_cores} /work/data/genomes/fly/Drosophila_melanogaster/UCSC/dm6/Sequence/BWAIndex/genome.fa ../${FILEID}.fastq > ${FILEID}.sai
bwa samse /work/data/genomes/fly/Drosophila_melanogaster/UCSC/dm6/Sequence/BWAIndex/genome.fa ${FILEID}.sai ../${FILEID}.fastq | samtools view -bS -o ${FILEID}.bam
samtools sort -m 2G -@ ${no_cores} -o ${FILEID}.s.bam ${FILEID}.bam
echo run time ${method} is $(expr `date +%s` - $start_time) s
mv ${FILEID}.s.bam ${FILEID}.${method}.bam
samtools index ${FILEID}.${method}.bam
rm ${FILEID}.bam
rm ${FILEID}.sai

### STAR aln
method="STAR"
start_time=`date +%s`
STAR --runThreadN 16 \
		--genomeDir /work/data/genomes/fly/Drosophila_melanogaster/UCSC/dm6/Sequence/STARgenomeIndex \
			--readFilesIn ../${FILEID}.fastq \
					--alignIntronMax 1 \
						--alignEndsType EndToEnd \
							--outSAMtype BAM SortedByCoordinate
echo run time ${method} is $(expr `date +%s` - $start_time) s
mv Aligned.sortedByCoord.out.bam ${FILEID}.${method}.bam
samtools index ${FILEID}.${method}.bam

### STAR aln
method="STAR_unique"
start_time=`date +%s`
STAR --runThreadN 16 \
		--genomeDir /work/data/genomes/fly/Drosophila_melanogaster/UCSC/dm6/Sequence/STARgenomeIndex \
			--readFilesIn ../${FILEID}.fastq \
				--outFilterMultimapNmax 1 \
					--alignIntronMax 1 \
						--alignEndsType EndToEnd \
							--outSAMtype BAM SortedByCoordinate
echo run time ${method} is $(expr `date +%s` - $start_time) s
mv Aligned.sortedByCoord.out.bam ${FILEID}.${method}.bam
samtools index ${FILEID}.${method}.bam

rm *.out*

#rm *.err; rm *.out; rm -rf *_out
