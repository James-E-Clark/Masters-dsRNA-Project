#!/bin/bash

ALN_DIR=STAR_Hilz_alignments
BED_DIR=Hilz_STAR_bed_files
CHR_BED=gencode.vM20.chromosomes.bed

mkdir -p ${BED_DIR}

# Rename samples
mv ${ALN_DIR}/SRR3659146.trimmedAligned.sortedByCoord.out.bam ${ALN_DIR}/day1_whole_testis_n1.sorted.bam
mv ${ALN_DIR}/SRR3659147.trimmedAligned.sortedByCoord.out.bam ${ALN_DIR}/day1_whole_testis_n2.sorted.bam
mv ${ALN_DIR}/SRR3659148.trimmedAligned.sortedByCoord.out.bam ${ALN_DIR}/day3_whole_testis_n1.sorted.bam
mv ${ALN_DIR}/SRR3659149.trimmedAligned.sortedByCoord.out.bam ${ALN_DIR}/day3_whole_testis_n2.sorted.bam
mv ${ALN_DIR}/SRR3659150.trimmedAligned.sortedByCoord.out.bam ${ALN_DIR}/day7_whole_testis_n1.sorted.bam
mv ${ALN_DIR}/SRR3659151.trimmedAligned.sortedByCoord.out.bam ${ALN_DIR}/day7_whole_testis_n2.sorted.bam

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

for i in day1 day3 day7
do

echo "Merging ${i} bed files..."
bedtools unionbedg -i ${BED_DIR}/${i}_whole_testis_n1.primary.bed ${BED_DIR}/${i}_whole_testis_n2.primary.bed > ${i}.tmp
# Sum the coverage columns
awk 'BEGIN{FS=OFS="\t"} {print $1, $2, $3, $4+$5}' ${i}.tmp >  ${BED_DIR}/${i}_whole_testis.merged.primary.bed
rm ${i}.tmp

done

# Merge all samples

bedtools unionbedg -i ${BED_DIR}/day1_whole_testis.merged.primary.bed ${BED_DIR}/day3_whole_testis.merged.primary.bed \
                      ${BED_DIR}/day7_whole_testis.merged.primary.bed > Hilz_all.merged.primary.tmp
awk 'BEGIN{FS=OFS="\t"} {print $1, $2, $3, $4+$5+$6}' Hilz_all.merged.primary.tmp >  ${BED_DIR}/Hilz_all.merged.primary.bed
rm Hilz_all.merged.primary.tmp
