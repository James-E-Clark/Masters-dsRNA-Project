#!/bin/bash

#SBATCH --mail-user=john.casement@ncl.ac.uk
#SBATCH --mail-type=FAIL,TIME_LIMIT
#SBATCH -n 11
#SBATCH -o logs/%x.%j.out
#SBATCH -A scbsu

module add Bowtie2
module add SAMtools

BASE_DIR=/nobackup/proj/scbsu/James_Clark
BT_IDX=${BASE_DIR}/genome/Hilz_bowtie_index
OUT_DIR=${BASE_DIR}/bowtie_Hilz_alignments

mkdir -p ${OUT_DIR}

bowtie2-build ${BASE_DIR}/genome/GRCm38.p6.genome.fa ${BASE_DIR}/genome/Hilz_bowtie_index

for i in ${BASE_DIR}/trimmed_Hilz_fastq/*.fastq
do

SAMPLE=$(basename ${i} .fastq)

# Run bowtie2 with settings for small RNA alignment
bowtie2 -p 11 -x ${BASE_DIR}/genome/Hilz_bowtie_index \
        -D 30 -R 3 -N 0 -L 12 -i S,1,0.50 \
        -U ${i} -S ${OUT_DIR}/${SAMPLE}.sam

# Convert to BAM, sort and index
samtools view -b ${OUT_DIR}/${SAMPLE}.sam -o ${OUT_DIR}/${SAMPLE}.bam
samtools sort ${OUT_DIR}/${SAMPLE}.bam -o ${OUT_DIR}/${SAMPLE}.sorted.bam
samtools index ${OUT_DIR}/${SAMPLE}.sorted.bam
rm ${OUT_DIR}/${SAMPLE}.sam

done 

