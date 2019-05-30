#!/bin/bash

mkdir -p HeLa_bed_files

# Rename samples
mv bowtie_HeLa_alignments/SRR5260520.sorted.bam bowtie_HeLa_alignments/Untreated.sorted.bam
mv bowtie_HeLa_alignments/SRR5260524.sorted.bam bowtie_HeLa_alignments/siCntrl.sorted.bam
mv bowtie_HeLa_alignments/SRR5260521.sorted.bam bowtie_HeLa_alignments/siSUV3_rep1.sorted.bam
mv bowtie_HeLa_alignments/SRR5260525.sorted.bam bowtie_HeLa_alignments/siSUV3_rep2.sorted.bam
mv bowtie_HeLa_alignments/SRR5260522.sorted.bam bowtie_HeLa_alignments/siPNPase_rep1.sorted.bam
mv bowtie_HeLa_alignments/SRR5260526.sorted.bam bowtie_HeLa_alignments/siPNPase_rep2.sorted.bam

for i in bowtie_HeLa_alignments/*.sorted.bam
do

SAMPLE=$(basename ${i} .sorted.bam)

# Remove any multi-mappers (XS field only present if more than one alignment was found for the read)
samtools view -H bowtie_HeLa_alignments/${SAMPLE}.sorted.bam > ${SAMPLE}.header.sam
samtools view -F 4 bowtie_HeLa_alignments/${SAMPLE}.sorted.bam | grep -v "XS:" | cat ${SAMPLE}.header.sam - | \
samtools view -b - > bowtie_HeLa_alignments/${SAMPLE}.unique.bam
samtools index bowtie_HeLa_alignments/${SAMPLE}.unique.bam
rm ${SAMPLE}.header.sam
rm bowtie_HeLa_alignments/${SAMPLE}.sorted.bam

echo "Running bedtools genomecov on sample ${SAMPLE}..."
# Get read coverage
bedtools genomecov -bg -ibam bowtie_HeLa_alignments/${SAMPLE}.unique.bam > HeLa_bed_files/${SAMPLE}.bed

echo "Running bedtools multicov on sample ${SAMPLE}..."
# Get reads per chromosome
bedtools multicov -bams bowtie_HeLa_alignments/${SAMPLE}.unique.bam -bed genome/gencode.v29.chromosomes.bed > HeLa_bed_files/${SAMPLE}.unique.reads.per.chr.bed

done

# Merge replicates

for i in siSUV3 siPNPase
do

echo "Merging ${i} bed files..."
bedtools unionbedg -i HeLa_bed_files/${i}_rep1.bed HeLa_bed_files/${i}_rep2.bed > ${i}.tmp
# Sum the coverage columns
awk 'BEGIN{FS=OFS="\t"} {print $1, $2, $3, $4+$5}' ${i}.tmp > HeLa_bed_files/${i}.merged.bed
rm ${i}.tmp

done

