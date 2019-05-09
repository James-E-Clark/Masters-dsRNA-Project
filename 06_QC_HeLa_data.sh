#!/bin/bash

mkdir -p HeLa_QC/FastQC

fastqc -o HeLa_QC/FastQC SRA_runs_HeLa/*.fastq
multiqc -o HeLa_QC/ HeLa_QC/FastQC/

