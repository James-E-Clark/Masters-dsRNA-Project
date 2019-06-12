#!/bin/bash

OUT_DIR=Hilz_intersections

mkdir ${OUT_DIR}

echo "Getting intersections between mouse F and Hilz..."
bedtools intersect -a mouse_STAR_bed_files/F_all.merged.bed -b Hilz_Bowtie2_bed_files/Hilz_all.merged.bed > Hilz_F_intersect.tmp
grep '^chr' Hilz_F_intersect.tmp > ${OUT_DIR}/Hilz_F_intersect.bed
bedtools intersect -wb -a ${OUT_DIR}/Hilz_F_intersect.bed -b genome/gencode.vM20.annotation.gff3 > ${OUT_DIR}/Hilz_F_intersect.annot.bed
rm Hilz_F_intersect.tmp

echo "Getting intersections between mouse R and Hilz..."
bedtools intersect -a mouse_STAR_bed_files/R_all.merged.bed -b Hilz_Bowtie2_bed_files/Hilz_all.merged.bed > Hilz_R_intersect.tmp
grep '^chr' Hilz_R_intersect.tmp > ${OUT_DIR}/Hilz_R_intersect.bed
bedtools intersect -wb -a ${OUT_DIR}/Hilz_R_intersect.bed -b genome/gencode.vM20.annotation.gff3 > ${OUT_DIR}/Hilz_R_intersect.annot.bed
rm Hilz_R_intersect.tmp
