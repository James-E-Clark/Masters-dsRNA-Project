#!/bin/bash

ALN_DIR=mouse_against_zebrafish_STAR_alignments
BED_DIR=mouse_against_zebrafish_STAR_bed_files

mkdir -p ${BED_DIR}

for i in ${ALN_DIR}/*Aligned.sortedByCoord.out.bam
do

SAMPLE=$(basename ${i} Aligned.sortedByCoord.out.bam)

# Remove any multi-mappers (keep reads denoted by NH:i:1 in STAR output)
samtools view -H ${ALN_DIR}/${SAMPLE}Aligned.sortedByCoord.out.bam > ${SAMPLE}.header.sam
samtools view -F 4 ${ALN_DIR}/${SAMPLE}Aligned.sortedByCoord.out.bam |  grep -w 'NH:i:1' | cat ${SAMPLE}.header.sam - | \
samtools view -b - > ${ALN_DIR}/${SAMPLE}.unique.bam
samtools index ${ALN_DIR}/${SAMPLE}.unique.bam
rm ${SAMPLE}.header.sam
#rm ${ALN_DIR}/${SAMPLE}Aligned.sortedByCoord.out.bam

echo "Running bedtools genomecov on sample ${SAMPLE}..."
# Get read coverage
bedtools genomecov -bg -ibam ${ALN_DIR}/${SAMPLE}.unique.bam > ${BED_DIR}/${SAMPLE}.bed

done

# Merge replicates

# Dicer KO samples are 2123, 2104 and 2464
# WT samples are 2501, 2128 and 2095
# Adult WT samples are A1 and A2

for i in F R
do

echo "Merging ${i} Dicer KO sample bed files..."
bedtools unionbedg -i ${BED_DIR}/${i}2123.bed ${BED_DIR}/${i}2104.bed ${BED_DIR}/${i}2464.bed > ${i}_Dicer_KO.tmp
# Sum the coverage columns
awk 'BEGIN{FS=OFS="\t"} {print $1, $2, $3, $4+$5+$6}' ${i}_Dicer_KO.tmp > ${BED_DIR}/${i}_Dicer_KO.merged.bed
rm ${i}_Dicer_KO.tmp

echo "Merging ${i} WT sample bed files..."
bedtools unionbedg -i ${BED_DIR}/${i}2501.bed ${BED_DIR}/${i}2128.bed ${BED_DIR}/${i}2095.bed > ${i}_WT.tmp
awk 'BEGIN{FS=OFS="\t"} {print $1, $2, $3, $4+$5+$6}' ${i}_WT.tmp > ${BED_DIR}/${i}_WT.merged.bed
rm ${i}_WT.tmp

echo "Merging ${i} Adult WT sample bed files..."
bedtools unionbedg -i ${BED_DIR}/${i}A1.bed ${BED_DIR}/${i}A2.bed > ${i}_Adult_WT.tmp
awk 'BEGIN{FS=OFS="\t"} {print $1, $2, $3, $4+$5}' ${i}_Adult_WT.tmp > ${BED_DIR}/${i}_Adult_WT.merged.bed
rm ${i}_Adult_WT.tmp

done
