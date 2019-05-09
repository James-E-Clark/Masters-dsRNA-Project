#!/bin/bash

mkdir -p bed_files

#for i in STAR_mm10_alignments/*.bam
for i in STAR_HeLa_Cell_alignments/*.bam
do

SAMPLE=$(basename ${i} Aligned.sortedByCoord.out.bam)

echo "Running bedtools genomecov on sample ${SAMPLE}..."
# output in bedgraph format
bedtools genomecov -bg -ibam ${i} > bed_files/${SAMPLE}.bedgraph

done
