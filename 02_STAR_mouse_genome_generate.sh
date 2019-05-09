#!/bin/bash

#SBATCH --mem-per-cpu=96000
#SBATCH -p bigmem
#SBATCH --mail-user=john.casement@ncl.ac.uk
#SBATCH --mail-type=FAIL,TIME_LIMIT
#SBATCH -n 11
#SBATCH -o logs/%x.%j.out
#SBATCH -A scbsu

STAR=/nobackup/proj/scbsu/software/STAR-2.6.0a/bin/Linux_x86_64/STAR
BASE_DIR=/nobackup/proj/scbsu/James_Clark
mkdir -p ${BASE_DIR}/genome/STAR_mm10_index

${STAR} --runThreadN 11 \
        --runMode genomeGenerate \
        --genomeDir ${BASE_DIR}/genome/STAR_mm10_index \
        --genomeFastaFiles ${BASE_DIR}/genome/GRCm38.p6.genome.fa \
        --sjdbGTFfile ${BASE_DIR}/genome/gencode.vM20.annotation.gtf \
        --sjdbOverhang 100 \
--limitGenomeGenerateRAM=96000000000
