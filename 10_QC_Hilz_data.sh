#!/bin/bash

mkdir -p Hilz_QC/FastQC

fastqc -o Hilz_QC/FastQC SRA_runs_Hilz/*.fastq
multiqc -o Hilz_QC/ Hilz_QC/FastQC/

mv Hilz_QC/multiqc_data  Hilz_QC/multiqc_data_raw/
mv Hilz_QC/multiqc_report.html Hilz_QC/multiqc_report_raw.html
