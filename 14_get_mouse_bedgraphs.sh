#!/bin/bash

mkdir -p mouse_bed_files

for i in STAR_mm10_alignments/*Aligned.sortedByCoord.out.bam
do

SAMPLE=$(basename ${i} Aligned.sortedByCoord.out.bam)

# Remove any multi-mappers (XS field only present if more than one alignment was found for the read)
samtools view -H STAR_mm10_alignments/${SAMPLE}Aligned.sortedByCoord.out.bam > ${SAMPLE}.header.sam
samtools view -F 4 STAR_mm10_alignments/${SAMPLE}Aligned.sortedByCoord.out.bam | grep -v "XS:" | cat ${SAMPLE}.header.sam - | \
samtools view -b - > STAR_mm10_alignments/${SAMPLE}.unique.bam
samtools index STAR_mm10_alignments/${SAMPLE}.unique.bam
rm ${SAMPLE}.header.sam
rm STAR_mm10_alignments/${SAMPLE}Aligned.sortedByCoord.out.bam

echo "Running bedtools genomecov on sample ${SAMPLE}..."
# Get read coverage
bedtools genomecov -bg -ibam STAR_mm10_alignments/${SAMPLE}.unique.bam > mouse_bed_files/${SAMPLE}.bed

echo "Running bedtools multicov on sample ${SAMPLE}..."
# Get reads per chromosome
bedtools multicov -bams STAR_mm10_alignments/${SAMPLE}.unique.bam -bed genome/gencode.vM20.chromosomes.bed > mouse_bed_files/${SAMPLE}.unique.reads.per.chr.bed

done

# Merge replicates

# Dicer KO samples are 2123, 2104 and 2464
# WT samples are 2501, 2128 and 2095
# Adult WT samples are A1 and A2

for i in F R
do

echo "Merging ${i} Dicer KO sample bed files..."
bedtools unionbedg -i mouse_bed_files/${i}2123.bed mouse_bed_files/${i}2104.bed mouse_bed_files/${i}2464.bed > ${i}_Dicer_KO.tmp
# Sum the coverage columns
awk 'BEGIN{FS=OFS="\t"} {print $1, $2, $3, $4+$5+$6}' ${i}_Dicer_KO.tmp > mouse_bed_files/${i}_Dicer_KO.merged.bed
rm ${i}_Dicer_KO.tmp

echo "Merging ${i} WT sample bed files..."
bedtools unionbedg -i mouse_bed_files/${i}2501.bed mouse_bed_files/${i}2128.bed mouse_bed_files/${i}2095.bed > ${i}_WT.tmp
awk 'BEGIN{FS=OFS="\t"} {print $1, $2, $3, $4+$5+$6}' ${i}_WT.tmp > mouse_bed_files/${i}_WT.merged.bed
rm ${i}_WT.tmp

echo "Merging ${i} Adult WT sample bed files..."
bedtools unionbedg -i mouse_bed_files/${i}A1.bed mouse_bed_files/${i}A2.bed > ${i}_Adult_WT.tmp
awk 'BEGIN{FS=OFS="\t"} {print $1, $2, $3, $4+$5}' ${i}_Adult_WT.tmp > mouse_bed_files/${i}_Adult_WT.merged.bed
rm ${i}_Adult_WT.tmp

done
