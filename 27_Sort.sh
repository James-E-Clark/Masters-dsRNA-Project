# Once SRA files have been aligned with STAR, the resulting BAM files are sorted:

module add SAMtools


for i in ../BAM/*.out.bam
do

SAMPLE=$(basename ${i} .bam)

samtools sort ../BAM/${SAMPLE}.bam -o ./${SAMPLE}.sorted.bam

echo "sample is ${SAMPLE}"
done
