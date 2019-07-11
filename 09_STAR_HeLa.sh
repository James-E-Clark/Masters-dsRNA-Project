#!/bin/bash

#SBATCH --mem-per-cpu=96000
#SBATCH -p bigmem
#SBATCH --mail-user=john.casement@ncl.ac.uk
#SBATCH --mail-type=FAIL,TIME_LIMIT
#SBATCH -n 11
#SBATCH -o logs/%x.%j.out
#SBATCH -A scbsu

module add SAMtools

STAR=/nobackup/proj/scbsu/software/STAR-2.6.0a/bin/Linux_x86_64/STAR

BASE_DIR=/nobackup/proj/scbsu/James_Clark
GENOME_DIR=/nobackup/proj/scbsu/James_Clark/genome/STAR_HeLa_index
OUT_DIR=${BASE_DIR}/STAR_HeLa_alignments

mkdir -p ${OUT_DIR}

for i in ${BASE_DIR}/SRA_runs_HeLa/*_1.fastq.gz
do

SAMPLE=$(basename ${i} _1.fastq.gz)

${STAR} --runThreadN 11 \
        --genomeDir ${GENOME_DIR} \
        --genomeLoad NoSharedMemory \
        --readFilesCommand zcat \
        --readFilesIn ${BASE_DIR}/SRA_runs_HeLa/${SAMPLE}_1.fastq.gz ${BASE_DIR}/SRA_runs_HeLa/${SAMPLE}_2.fastq.gz \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix ${OUT_DIR}/${SAMPLE} \
        --twopassMode Basic \
        --twopass1readsN -1

# Create index file
samtools index ${OUT_DIR}/${SAMPLE}Aligned.sortedByCoord.out.bam

done

