#!/bin/bash

# Coverage bedgraphs for Hilz data (https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP076460)

BAM_DIR=/home/john/customers/Werner/James_Clark/bowtie_Hilz_alignments
BG_DIR=/home/john/customers/Werner/James_Clark/bed_files

mkdir -p bed_files

# Rename files
mv ${BAM_DIR}/SRR3659146.trimmed.sorted.bam ${BAM_DIR}/day1_whole_testis_n1.bam
mv ${BAM_DIR}/SRR3659147.trimmed.sorted.bam ${BAM_DIR}/day1_whole_testis_n2.bam
mv ${BAM_DIR}/SRR3659148.trimmed.sorted.bam ${BAM_DIR}/day3_whole_testis_n1.bam
mv ${BAM_DIR}/SRR3659149.trimmed.sorted.bam ${BAM_DIR}/day3_whole_testis_n2.bam
mv ${BAM_DIR}/SRR3659150.trimmed.sorted.bam ${BAM_DIR}/day7_whole_testis_n1.bam
mv ${BAM_DIR}/SRR3659151.trimmed.sorted.bam ${BAM_DIR}/day7_whole_testis_n2.bam

for i in ${BAM_DIR}/*.bam
do

SAMPLE=$(basename ${i} .bam)

echo "Running bedtools genomecov on sample ${SAMPLE}..."
# output in bedgraph format
bedtools genomecov -bg -ibam ${i} > ${BG_DIR}/${SAMPLE}.bedgraph

done

# Merge replicates
echo "Merging day 1 bedgraphs..."
bedtools unionbedg -i ${BG_DIR}/day1_whole_testis_n1.bedgraph ${BG_DIR}/day1_whole_testis_n2.bedgraph > day1.tmp
# Sum the coverage columns
awk 'BEGIN{FS=OFS="\t"} {print $1, $2, $3, $4+$5}' day1.tmp > ${BG_DIR}/day1_whole_testis.union.bedgraph
rm day1.tmp

echo "Merging day 3 bedgraphs..."
bedtools unionbedg -i ${BG_DIR}/day3_whole_testis_n1.bedgraph ${BG_DIR}/day3_whole_testis_n2.bedgraph > day3.tmp
awk 'BEGIN{FS=OFS="\t"} {print $1, $2, $3, $4+$5}' day3.tmp > ${BG_DIR}/day3_whole_testis.union.bedgraph
rm day3.tmp

echo "Merging day 7 bedgraphs..."
bedtools unionbedg -i ${BG_DIR}/day7_whole_testis_n1.bedgraph ${BG_DIR}/day7_whole_testis_n2.bedgraph > day7.tmp
awk 'BEGIN{FS=OFS="\t"} {print $1, $2, $3, $4+$5}' day7.tmp > ${BG_DIR}/day7_whole_testis.union.bedgraph
rm day7.tmp

