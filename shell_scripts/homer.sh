#!/bin/sh
# homer.sh

# first argument sample table
TABLE="$1" 

# collect file names
FILES=`cut -f 1 $TABLE `

# make tag directory for bam files
for FILE in $FILES
do
	echo $FILE
	makeTagDirectory ${FILE}.dir ${FILE}.bam
done

# iterate through samples
SAMPLES=`cut -f 2 $TABLE | uniq`

for SAMPLE in $SAMPLES
do
  # get input
	INPUT=`cat $TABLE | grep $SAMPLE | grep "INPUT" | cut -f 1`

  # get non-inputs
	ASSAYS=`cat $TABLE | grep $SAMPLE | grep -v "INPUT" | cut -f 1`
  
	for ASSAY in $ASSAYS
	do
	  # show ChIP - Input pairs
		echo $ASSAY $INPUT

    # create total reads and Input normalized bedgraph
		makeUCSCfile ${ASSAY}.dir/  -i ${INPUT}.dir/  -o ${ASSAY}.INPnorm.bedgraph
  
    # peak finding against input
		findPeaks ${ASSAY}.dir/  -i ${INPUT}.dir/ -style factor -F 2 -o ${ASSAY}.factor.F2.txt
    
    # convert txt to BED
		cat ${ASSAY}.factor.F2.txt | grep -v "#" | cut -f 2,3,4 > ${ASSAY}.factor.F2.bed
	done

done