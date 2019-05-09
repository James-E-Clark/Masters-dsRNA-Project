#!/bin/bash

# Sample relations from SRA:
# SRR5260520	Untreated_J2_IP
# SRR5260521	siSUV3_J2_IP (rep1)
# SRR5260522	siPNPase_J2_IP (rep1)
# SRR5260524	siCntrl_J2_IP
# SRR5260525	siSUV3_J2_IP (rep2)
# SRR5260526	siPNPase_J2_IP (rep2)

BG_DIR=/home/john/customers/Werner/James_Clark/bed_files

# Rename bedgraphs
mv ${BG_DIR}/SRR5260520.bedgraph ${BG_DIR}/Untreated.bedgraph
mv ${BG_DIR}/SRR5260524.bedgraph ${BG_DIR}/siCntrl.bedgraph
mv ${BG_DIR}/SRR5260521.bedgraph ${BG_DIR}/siSUV3_rep1.bedgraph
mv ${BG_DIR}/SRR5260525.bedgraph ${BG_DIR}/siSUV3_rep2.bedgraph
mv ${BG_DIR}/SRR5260522.bedgraph ${BG_DIR}/siPNPase_rep1.bedgraph
mv ${BG_DIR}/SRR5260526.bedgraph ${BG_DIR}/siPNPase_rep2.bedgraph

# Merge by sample type
echo "Merging Untreated and siCntrl bedgraphs..."
bedtools unionbedg -i ${BG_DIR}/Untreated.bedgraph ${BG_DIR}/siCntrl.bedgraph > Cntrl_Untreated.tmp
# Sum the coverage columns
awk 'BEGIN{FS=OFS="\t"} {print $1, $2, $3, $4+$5}' Cntrl_Untreated.tmp > ${BG_DIR}/siCntrl_Untreated.union.bedgraph
rm Cntrl_Untreated.tmp

echo "Merging siSUV3 bedgraphs..."
bedtools unionbedg -i ${BG_DIR}/siSUV3_rep1.bedgraph ${BG_DIR}/siSUV3_rep2.bedgraph > siSUV3.tmp
# Sum the coverage columns
awk 'BEGIN{FS=OFS="\t"} {print $1, $2, $3, $4+$5}' siSUV3.tmp > ${BG_DIR}/siSUV3.union.bedgraph
rm siSUV3.tmp

echo "Merging siPNPase bedgraphs..."
bedtools unionbedg -i ${BG_DIR}/siPNPase_rep1.bedgraph ${BG_DIR}/siPNPase_rep2.bedgraph > siPNPase.tmp
# Sum the coverage columns
awk 'BEGIN{FS=OFS="\t"} {print $1, $2, $3, $4+$5}' siPNPase.tmp > ${BG_DIR}/siPNPase.union.bedgraph
rm siPNPase.tmp

echo "Done"
