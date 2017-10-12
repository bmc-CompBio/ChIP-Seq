#! /bin/bash

module load ngs/sratoolkit/2.8.0

prefetch SRR495368 #IP
prefetch SRR495378 #input

fastq-dump SRR495368
fastq-dump SRR495378

FILES=($(ls -1 *.fastq))
# get size of array
NUMFASTQ=${#FILES[@]}
 
# now submit to SLURM
if [ $NUMFASTQ -ge 0 ]; then
	sbatch --array=1-$NUMFASTQ alignTest_arrays.sbatch
fi

while read method junk
do
	sbatch --error=homer_${method}.err --output=homer_${method}.out --export=method=${method} homer.sbatch 
done < methods.txt
