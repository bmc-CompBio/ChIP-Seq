#! /bin/bash
#
# homer.sbatch
#
#SBATCH -J homer                     # A single job name for the array
#SBATCH -p slim18                    # Partition
#SBATCH -n 9                         # 18 cores
#SBATCH -N 1                         # nodes required
#SBATCH -t 0-01:00                   # running time
#SBATCH --mem 10000                  # Memory request
#SBATCH -o homer_%A_%a.out       # Standard output
#SBATCH -e homer_%A_%a.err       # Standard error

module load ngs/Homer/4.9
module load ngs/samtools/1.3.1

makeTagDirectory ${method}_tag/ip/ SRR495368_out/SRR495368.${method}.bam -single
makeTagDirectory ${method}_tag/input/ SRR495378_out/SRR495378.${method}.bam -single

findPeaks ${method}_tag/ip/ -i ${method}_tag/input/ \
	-style factor -fragLength 200 -inputFragLength 200 \
	-F 2 -C 0 > ${method}_peaks.txt
pos2bed.pl ${method}_peaks.txt > ${method}_peaks.bed

rm -rf ${method}_tag
