#!/bin/bash

# Short reads with many small RNA adaptors, so trim these out with trimmomatic

#SBATCH --mem-per-cpu=5000
#SBATCH --mail-user=john.casement@ncl.ac.uk
#SBATCH --mail-type=FAIL,TIME_LIMIT
#SBATCH -n 11
#SBATCH -o logs/%x.%j.out
#SBATCH -A scbsu

module load Java/1.8.0_144

BASE_DIR=/nobackup/proj/scbsu/James_Clark
OUT_DIR=${BASE_DIR}/trimmed_Hilz_fastq
mkdir -p ${OUT_DIR}

for i in ${BASE_DIR}/SRA_runs_Hilz/*.fastq
do

SAMPLE=$(basename ${i} .fastq)

# Check if trimmed file already exists 
if [ ! -f ${OUT_DIR}/${SAMPLE}.trimmed.fastq ]
then

echo "Trimming ${SAMPLE}.fastq"
java -jar /nobackup/proj/scbsu/software/Trimmomatic/Trimmomatic-0.36/trimmomatic-0.36.jar \
          SE \
         -threads 11 \
         -phred33 \
         -trimlog ${OUT_DIR}/${SAMPLE}.trimlog.txt \
         ${BASE_DIR}/SRA_runs_Hilz/${SAMPLE}.fastq \
         ${OUT_DIR}/${SAMPLE}.trimmed.fastq \
         ILLUMINACLIP:/nobackup/proj/scbsu/software/Trimmomatic/Trimmomatic-0.36/adapters/all_adapters.fa:1:5:5:8:true \
         AVGQUAL:20 \
         MINLEN:15

fi

done
 
