#!/bin/bash

ALN_DIR=HeLa_Bowtie2_alignments
BED_DIR=HeLa_Bowtie2_bed_files
CHR_BED=gencode.v29.chromosomes.bed

mkdir -p ${BED_DIR}

# Rename samples
mv ${ALN_DIR}/SRR5260520.sorted.bam ${ALN_DIR}/Untreated.sorted.bam
mv ${ALN_DIR}/SRR5260524.sorted.bam ${ALN_DIR}/siCntrl.sorted.bam
mv ${ALN_DIR}/SRR5260521.sorted.bam ${ALN_DIR}/siSUV3_rep1.sorted.bam
mv ${ALN_DIR}/SRR5260525.sorted.bam ${ALN_DIR}/siSUV3_rep2.sorted.bam
mv ${ALN_DIR}/SRR5260522.sorted.bam ${ALN_DIR}/siPNPase_rep1.sorted.bam
mv ${ALN_DIR}/SRR5260526.sorted.bam ${ALN_DIR}/siPNPase_rep2.sorted.bam

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

for i in siSUV3 siPNPase
do

echo "Merging ${i} bed files..."
bedtools unionbedg -i ${BED_DIR}/${i}_rep1.bed ${BED_DIR}/${i}_rep2.bed > ${i}.tmp
# Sum the coverage columns
awk 'BEGIN{FS=OFS="\t"} {print $1, $2, $3, $4+$5}' ${i}.tmp > ${BED_DIR}/${i}.merged.bed
rm ${i}.tmp

done
