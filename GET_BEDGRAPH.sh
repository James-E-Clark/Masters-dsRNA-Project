#Once the BAM files are sorted, a bedgraph can be created

#ROCKET COMMANDS:
#!/bin/bash

#SBATCH --mem-per-cpu=5000
#SBATCH --mail-user=jamesclark1996@gmail.com
#SBATCH --mail-type=FAIL,TIME_LIMIT
#SBATCH -n 11
#SBATCH -o logs/%x.%j.out
#SBATCH -A awed

#ADD BEDTOOLS PACKAGE:
module add BEDTools

#DEFINE FILEPATHS:
ALN_DIR=/nobackup/proj/awed/james_clark/Gao_et_al/BEDGRAPHS
BED_DIR=/nobackup/proj/awed/james_clark/Gao_et_al/BEDGRAPHS/BEDGRAPH
CHR_BED=gencode.vM20.chromosomes.bed

mkdir -p ${BED_DIR}

#FOR LOOP TO DO ALL FILES CONSECUTIVELY:
for i in ${ALN_DIR}/*.sorted.bam
do

SAMPLE=$(basename ${i} .Aligned.sortedByCoord.out.sorted.bam)
bedtools genomecov -bg -ibam ${ALN_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.sorted.bam > ${BED_DIR}/${SAMPLE}.bed
#echo "sample is ${SAMPLE}"
done

#FINISHED. (The hashed out line above #echo "sample is ${SAMPLE}" can be used to check that the sample name is what is expected, unhash the line and run the scipt to have the sample names printed on the screen. 
