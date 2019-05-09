#!/bin/bash

FASTQ_DUMP=sratoolkit.2.9.6-ubuntu64/bin/fastq-dump

for i in SRR3659146 SRR3659147 SRR3659148 SRR3659149 SRR3659150 SRR3659151
do

${FASTQ_DUMP} --outdir SRA_runs_Hilz ${i}

done

