#!/bin/bash

#SBATCH --mail-user=john.casement@ncl.ac.uk
#SBATCH --mail-type=FAIL,TIME_LIMIT
#SBATCH -n 11
#SBATCH -o logs/%x.%j.out
#SBATCH -A scbsu

module add Bowtie2
module add SAMtools

BASE_DIR=/nobackup/proj/scbsu/James_Clark
BT_IDX=${BASE_DIR}/genome/HeLa_bowtie_index
OUT_DIR=${BASE_DIR}/bowtie_HeLa_alignments

mkdir -p ${OUT_DIR}

bowtie2-build ${BASE_DIR}/genome/GRCh38.primary_assembly.genome.fa ${BASE_DIR}/genome/HeLa_bowtie_index

for i in ${BASE_DIR}/SRA_runs_HeLa/*_1.fastq.gz
do

SAMPLE=$(basename ${i} _1.fastq.gz)

bowtie2 -p 11 -x ${BASE_DIR}/genome/HeLa_bowtie_index --very-sensitive \
        -1 ${BASE_DIR}/SRA_runs_HeLa/${SAMPLE}_1.fastq.gz \
        -2 ${BASE_DIR}/SRA_runs_HeLa/${SAMPLE}_2.fastq.gz \
        -S ${OUT_DIR}/${SAMPLE}.sam

# Convert to BAM, sort and index
samtools view -b ${OUT_DIR}/${SAMPLE}.sam -o ${OUT_DIR}/${SAMPLE}.bam
samtools sort ${OUT_DIR}/${SAMPLE}.bam -o ${OUT_DIR}/${SAMPLE}.sorted.bam
samtools index ${OUT_DIR}/${SAMPLE}.sorted.bam
rm ${OUT_DIR}/${SAMPLE}.sam

done 

