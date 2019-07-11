#!/bin/bash

FASTQ_DUMP=sratoolkit.2.9.6-ubuntu64/bin/fastq-dump

for i in SRR5260520 SRR5260521 SRR5260522 SRR5260524 SRR5260525 SRR5260526
do

${FASTQ_DUMP} --gzip --split-files --outdir SRA_runs_HeLa ${i}

done

