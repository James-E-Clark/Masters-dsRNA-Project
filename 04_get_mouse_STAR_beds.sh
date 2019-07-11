#!/bin/bash

ALN_DIR=STAR_alignments
BED_DIR=mouse_STAR_bed_files
CHR_BED=gencode.vM20.chromosomes.bed

mkdir -p ${BED_DIR}

for i in ${ALN_DIR}/*Aligned.sortedByCoord.out.bam
do

SAMPLE=$(basename ${i} Aligned.sortedByCoord.out.bam)

## Remove any multi-mappers (keep only primary alignments with SAM flag 256)
samtools view -H ${ALN_DIR}/${SAMPLE}Aligned.sortedByCoord.out.bam > ${SAMPLE}.header.sam
samtools view -F 256 ${ALN_DIR}/${SAMPLE}Aligned.sortedByCoord.out.bam | cat ${SAMPLE}.header.sam - > ${SAMPLE}.primary.sam
samtools view -b ${SAMPLE}.primary.sam > ${ALN_DIR}/${SAMPLE}.primary.bam
samtools index ${ALN_DIR}/${SAMPLE}.primary.bam
rm ${SAMPLE}.header.sam
rm ${SAMPLE}.primary.sam

echo "Running bedtools genomecov on sample ${SAMPLE}..."
# Get read coverage
bedtools genomecov -bg -ibam ${ALN_DIR}/${SAMPLE}.primary.bam > ${BED_DIR}/${SAMPLE}.primary.bed

done

# Merge replicates

# Dicer KO samples are 2123, 2104 and 2464
# WT samples are 2501, 2128 and 2095
# Adult WT samples are A1 and A2

for i in F R
do

echo "Merging ${i} Dicer KO sample bed files..."
bedtools unionbedg -i ${BED_DIR}/${i}2123.primary.bed ${BED_DIR}/${i}2104.primary.bed ${BED_DIR}/${i}2464.primary.bed > ${i}_Dicer_KO.tmp
# Sum the coverage columns
awk 'BEGIN{FS=OFS="\t"} {print $1, $2, $3, $4+$5+$6}' ${i}_Dicer_KO.tmp > ${BED_DIR}/${i}_Dicer_KO.merged.primary.bed
rm ${i}_Dicer_KO.tmp

echo "Merging ${i} WT sample bed files..."
bedtools unionbedg -i ${BED_DIR}/${i}2501.primary.bed ${BED_DIR}/${i}2128.primary.bed ${BED_DIR}/${i}2095.primary.bed > ${i}_WT.tmp
awk 'BEGIN{FS=OFS="\t"} {print $1, $2, $3, $4+$5+$6}' ${i}_WT.tmp > ${BED_DIR}/${i}_WT.merged.primary.bed
rm ${i}_WT.tmp

echo "Merging ${i} Adult WT sample bed files..."
bedtools unionbedg -i ${BED_DIR}/${i}A1.primary.bed ${BED_DIR}/${i}A2.primary.bed > ${i}_Adult_WT.tmp
awk 'BEGIN{FS=OFS="\t"} {print $1, $2, $3, $4+$5}' ${i}_Adult_WT.tmp > ${BED_DIR}/${i}_Adult_WT.merged.primary.bed
rm ${i}_Adult_WT.tmp

echo "Merging ${i} all sample bed files..."
bedtools unionbedg -i ${BED_DIR}/${i}_Dicer_KO.merged.primary.bed ${BED_DIR}/${i}_WT.merged.primary.bed ${BED_DIR}/${i}_Adult_WT.merged.primary.bed  > ${i}_all.tmp
awk 'BEGIN{FS=OFS="\t"} {print $1, $2, $3, $4+$5+$6}' ${i}_all.tmp > ${i}_all.merged.primary.bed
rm ${i}_all.tmp

done
