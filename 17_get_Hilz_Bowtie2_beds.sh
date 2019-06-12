#!/bin/bash

ALN_DIR=Hilz_Bowtie2_alignments
BED_DIR=Hilz_Bowtie2_bed_files
CHR_BED=gencode.vM20.chromosomes.bed

mkdir -p ${BED_DIR}

# Rename samples
mv ${ALN_DIR}/SRR3659146.trimmed.sorted.bam ${ALN_DIR}/day1_whole_testis_n1.sorted.bam
mv ${ALN_DIR}/SRR3659147.trimmed.sorted.bam ${ALN_DIR}/day1_whole_testis_n2.sorted.bam
mv ${ALN_DIR}/SRR3659148.trimmed.sorted.bam ${ALN_DIR}/day3_whole_testis_n1.sorted.bam
mv ${ALN_DIR}/SRR3659149.trimmed.sorted.bam ${ALN_DIR}/day3_whole_testis_n2.sorted.bam
mv ${ALN_DIR}/SRR3659150.trimmed.sorted.bam ${ALN_DIR}/day7_whole_testis_n1.sorted.bam
mv ${ALN_DIR}/SRR3659151.trimmed.sorted.bam ${ALN_DIR}/day7_whole_testis_n2.sorted.bam


for i in ${ALN_DIR}/*.sorted.bam
do

SAMPLE=$(basename ${i} .sorted.bam)

# Remove any multi-mappers (XS field only present if more than one alignment was found for the read)
samtools view -H ${ALN_DIR}/${SAMPLE}.sorted.bam > ${SAMPLE}.header.sam
samtools view -F 4 ${ALN_DIR}/${SAMPLE}.sorted.bam | grep -v "XS:" | cat ${SAMPLE}.header.sam - | \
samtools view -b - > ${ALN_DIR}/${SAMPLE}.unique.bam
samtools index ${ALN_DIR}/${SAMPLE}.unique.bam
rm ${SAMPLE}.header.sam
#rm ${ALN_DIR}/${SAMPLE}.sorted.bam

echo "Running bedtools genomecov on sample ${SAMPLE}..."
# Get read coverage
bedtools genomecov -bg -ibam ${ALN_DIR}/${SAMPLE}.unique.bam > ${BED_DIR}/${SAMPLE}.bed

echo "Running bedtools multicov on sample ${SAMPLE}..."
# Get reads per chromosome
bedtools multicov -bams ${ALN_DIR}/${SAMPLE}.unique.bam -bed genome/${CHR_BED} > ${BED_DIR}/${SAMPLE}.unique.reads.per.chr.bed

done

# Merge replicates

for i in day1 day3 day7
do

echo "Merging ${i} bed files..."
bedtools unionbedg -i ${BED_DIR}/${i}_whole_testis_n1.bed ${BED_DIR}/${i}_whole_testis_n2.bed > ${i}.tmp
# Sum the coverage columns
awk 'BEGIN{FS=OFS="\t"} {print $1, $2, $3, $4+$5}' ${i}.tmp >  ${BED_DIR}/${i}_whole_testis.merged.bed
rm ${i}.tmp

done

# Merge all samples

bedtools unionbedg -i ${BED_DIR}/day1_whole_testis.merged.bed ${BED_DIR}/day3_whole_testis.merged.bed \
                      ${BED_DIR}/day7_whole_testis.merged.bed > Hilz_all.merged.tmp
awk 'BEGIN{FS=OFS="\t"} {print $1, $2, $3, $4+$5+$6}' Hilz_all.merged.tmp >  ${BED_DIR}/Hilz_all.merged.bed
rm Hilz_all.merged.tmp


