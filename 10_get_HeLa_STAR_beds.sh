#!/bin/bash

ALN_DIR=STAR_HeLa_alignments
BED_DIR=HeLa_STAR_bed_files
CHR_BED=gencode.v29.chromosomes.bed

mkdir -p ${BED_DIR}

# Rename samples
mv ${ALN_DIR}/SRR5260520Aligned.sortedByCoord.out.bam ${ALN_DIR}/Untreated.sorted.bam
mv ${ALN_DIR}/SRR5260524Aligned.sortedByCoord.out.bam ${ALN_DIR}/siCntrl.sorted.bam
mv ${ALN_DIR}/SRR5260521Aligned.sortedByCoord.out.bam ${ALN_DIR}/siSUV3_rep1.sorted.bam
mv ${ALN_DIR}/SRR5260525Aligned.sortedByCoord.out.bam ${ALN_DIR}/siSUV3_rep2.sorted.bam
mv ${ALN_DIR}/SRR5260522Aligned.sortedByCoord.out.bam ${ALN_DIR}/siPNPase_rep1.sorted.bam
mv ${ALN_DIR}/SRR5260526Aligned.sortedByCoord.out.bam ${ALN_DIR}/siPNPase_rep2.sorted.bam

for i in ${ALN_DIR}/*.sorted.bam
do

SAMPLE=$(basename ${i} .sorted.bam)

## Remove any multi-mappers (keep only primary alignments with SAM flag 256)
samtools view -H ${ALN_DIR}/${SAMPLE}.sorted.bam > ${SAMPLE}.header.sam
samtools view -F 256 ${ALN_DIR}/${SAMPLE}.sorted.bam | cat ${SAMPLE}.header.sam - > ${SAMPLE}.primary.sam
samtools view -b ${SAMPLE}.primary.sam > ${ALN_DIR}/${SAMPLE}.primary.bam
samtools index ${ALN_DIR}/${SAMPLE}.primary.bam
rm ${SAMPLE}.header.sam
rm ${SAMPLE}.primary.sam

echo "Running bedtools genomecov on sample ${SAMPLE}..."
# Get read coverage
bedtools genomecov -bg -ibam ${ALN_DIR}/${SAMPLE}.primary.bam > ${BED_DIR}/${SAMPLE}.primary.bed

done

# Merge replicates

for i in siSUV3 siPNPase
do

echo "Merging ${i} bed files..."
bedtools unionbedg -i ${BED_DIR}/${i}_rep1.primary.bed ${BED_DIR}/${i}_rep2.primary.bed > ${i}.tmp
# Sum the coverage columns
awk 'BEGIN{FS=OFS="\t"} {print $1, $2, $3, $4+$5}' ${i}.tmp > ${BED_DIR}/${i}.merged.primary.bed
rm ${i}.tmp

done

