#!/bin/bash

# Dicer KO samples are 2123, 2104 and 2464
# WT samples are 2501, 2128 and 2095
# Adult WT samples are A1 and A2

BG_DIR=/home/john/customers/Werner/James_Clark/bed_files

for i in F R
do

echo "Merging ${i} Dicer KO sample bedgraphs..."
bedtools unionbedg -i ${BG_DIR}/${i}2123.bedgraph ${BG_DIR}/${i}2104.bedgraph ${BG_DIR}/${i}2464.bedgraph > ${i}_Dicer_KO.tmp
# Sum the coverage columns
awk 'BEGIN{FS=OFS="\t"} {print $1, $2, $3, $4+$5+$6}' ${i}_Dicer_KO.tmp > ${BG_DIR}/${i}_Dicer_KO.union.bedgraph
rm ${i}_Dicer_KO.tmp

echo "Merging ${i} WT sample bedgraphs..."
bedtools unionbedg -i ${BG_DIR}/${i}2501.bedgraph ${BG_DIR}/${i}2128.bedgraph ${BG_DIR}/${i}2095.bedgraph > ${i}_WT.tmp
awk 'BEGIN{FS=OFS="\t"} {print $1, $2, $3, $4+$5+$6}' ${i}_WT.tmp > ${BG_DIR}/${i}_WT.union.bedgraph
rm ${i}_WT.tmp

echo "Merging ${i} Adult WT sample bedgraphs..."
bedtools unionbedg -i ${BG_DIR}/${i}A1.bedgraph ${BG_DIR}/${i}A2.bedgraph > ${i}_Adult_WT.tmp
# Sum the coverage columns
awk 'BEGIN{FS=OFS="\t"} {print $1, $2, $3, $4+$5}' ${i}_Adult_WT.tmp > ${BG_DIR}/${i}_Adult_WT.union.bedgraph
rm ${i}_Adult_WT.tmp

done

echo "Done"
