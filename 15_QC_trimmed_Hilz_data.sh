#!/bin/bash

mkdir -p Hilz_QC/trimmed_FastQC

fastqc -o Hilz_QC/trimmed_FastQC SRA_runs_Hilz/*.trimmed.fastq
multiqc -o Hilz_QC Hilz_QC/trimmed_FastQC/

mv Hilz_QC/multiqc_data  Hilz_QC/multiqc_data_trimmed
mv Hilz_QC/multiqc_report.html Hilz_QC/multiqc_report_trimmed.html

