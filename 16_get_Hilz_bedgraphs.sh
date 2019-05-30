#!/bin/bash

mkdir -p Hilz_bed_files

# Rename samples
mv bowtie_Hilz_alignments/SRR3659146.trimmed.sorted.bam bowtie_Hilz_alignments/day1_whole_testis_n1.sorted.bam
mv bowtie_Hilz_alignments/SRR3659147.trimmed.sorted.bam bowtie_Hilz_alignments/day1_whole_testis_n2.sorted.bam
mv bowtie_Hilz_alignments/SRR3659148.trimmed.sorted.bam bowtie_Hilz_alignments/day3_whole_testis_n1.sorted.bam
mv bowtie_Hilz_alignments/SRR3659149.trimmed.sorted.bam bowtie_Hilz_alignments/day3_whole_testis_n2.sorted.bam
mv bowtie_Hilz_alignments/SRR3659150.trimmed.sorted.bam bowtie_Hilz_alignments/day7_whole_testis_n1.sorted.bam
mv bowtie_Hilz_alignments/SRR3659151.trimmed.sorted.bam bowtie_Hilz_alignments/day7_whole_testis_n2.sorted.bam


for i in bowtie_Hilz_alignments/*.sorted.bam
do

SAMPLE=$(basename ${i} .sorted.bam)

# Remove any multi-mappers (XS field only present if more than one alignment was found for the read)
samtools view -H bowtie_Hilz_alignments/${SAMPLE}.sorted.bam > ${SAMPLE}.header.sam
samtools view -F 4 bowtie_Hilz_alignments/${SAMPLE}.sorted.bam | grep -v "XS:" | cat ${SAMPLE}.header.sam - | \
samtools view -b - > bowtie_Hilz_alignments/${SAMPLE}.unique.bam
samtools index bowtie_Hilz_alignments/${SAMPLE}.unique.bam
rm ${SAMPLE}.header.sam
rm bowtie_Hilz_alignments/${SAMPLE}.sorted.bam

echo "Running bedtools genomecov on sample ${SAMPLE}..."
# Get read coverage
bedtools genomecov -bg -ibam bowtie_Hilz_alignments/${SAMPLE}.unique.bam > Hilz_bed_files/${SAMPLE}.bed

echo "Running bedtools multicov on sample ${SAMPLE}..."
# Get reads per chromosome
bedtools multicov -bams bowtie_Hilz_alignments/${SAMPLE}.unique.bam -bed genome/gencode.vM20.chromosomes.bed > Hilz_bed_files/${SAMPLE}.unique.reads.per.chr.bed

done

# Merge replicates

for i in day1 day3 day7
do

echo "Merging ${i} bed files..."
bedtools unionbedg -i Hilz_bed_files/${i}_whole_testis_n1.bed Hilz_bed_files/${i}_whole_testis_n2.bed > ${i}.tmp
# Sum the coverage columns
awk 'BEGIN{FS=OFS="\t"} {print $1, $2, $3, $4+$5}' ${i}.tmp >  Hilz_bed_files/${i}_whole_testis.merged.bed
rm ${i}.tmp

done
