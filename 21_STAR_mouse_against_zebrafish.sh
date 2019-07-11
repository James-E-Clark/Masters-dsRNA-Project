#!/bin/bash

BASE_DIR=/home/john/customers/Werner
OUT_DIR=${BASE_DIR}/STAR_zebrafish_alignments

mkdir -p ${BASE_DIR}/zebrafish_genome/STAR_danrer_index
mkdir -p ${OUT_DIR}

# Generate genome
STAR --runThreadN 2 \
     --runMode genomeGenerate \
     --genomeDir ${BASE_DIR}/zebrafish_genome/STAR_danrer_index \
     --genomeFastaFiles ${BASE_DIR}/zebrafish_genome/zf-spike.fa

# Run STAR

for i in ${BASE_DIR}/Natasya/fastq/*R1.fq.gz.PwU.qtrim.fq.gz
do

SAMPLE=$(echo ${i} | awk -F "_R|/" '{print $(NF-1)}')

STAR --runThreadN 2 \
     --genomeDir ${BASE_DIR}/zebrafish_genome/STAR_danrer_index \
     --genomeLoad NoSharedMemory \
     --readFilesCommand zcat \
     --readFilesIn ${BASE_DIR}/Natasya/fastq/${SAMPLE}_R1.fq.gz.PwU.qtrim.fq.gz ${BASE_DIR}/Natasya/fastq/${SAMPLE}_R2.fq.gz.PwU.qtrim.fq.gz \
     --outSAMtype BAM SortedByCoordinate \
     --outFileNamePrefix ${OUT_DIR}/${SAMPLE} \
     --twopassMode Basic \
     --twopass1readsN -1

# Create index file
samtools index ${OUT_DIR}/${SAMPLE}Aligned.sortedByCoord.out.bam

done

