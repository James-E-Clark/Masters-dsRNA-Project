#!/bin/bash

#SBATCH --mem-per-cpu=96000
#SBATCH -p bigmem
#SBATCH --mail-user=john.casement@ncl.ac.uk
#SBATCH --mail-type=FAIL,TIME_LIMIT
#SBATCH -n 11
#SBATCH -o logs/%x.%j.out
#SBATCH -A scbsu

STAR=/nobackup/proj/scbsu/software/STAR-2.6.0a/bin/Linux_x86_64/STAR
module add SAMtools

BASE_DIR=/nobackup/proj/scbsu/James_Clark
GENOME_DIR=${BASE_DIR}/genome/STAR_Hilz_index
OUT_DIR=${BASE_DIR}/STAR_Hilz_alignments

mkdir ${OUT_DIR}

for i in ${BASE_DIR}/trimmed_Hilz_fastq/*.fastq
do

SAMPLE=$(basename ${i} .fastq)

${STAR} --runThreadN 11 \
        --genomeDir ${GENOME_DIR} \
        --genomeLoad NoSharedMemory \
        --readFilesIn ${BASE_DIR}/trimmed_Hilz_fastq/${SAMPLE}.fastq \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix ${OUT_DIR}/${SAMPLE} \
        --twopassMode Basic \
        --twopass1readsN -1

# Create index file
samtools index ${OUT_DIR}/${SAMPLE}Aligned.sortedByCoord.out.bam

done

