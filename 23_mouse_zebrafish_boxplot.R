# Boxplots for mouse data

setwd("~/customers/Werner/James_Clark/")

library(ggplot2)

# F samples
F_Dicer_KO = read.table("mouse_against_zebrafish_STAR_bed_files/F_Dicer_KO.merged.bed",
                        sep = "\t", header = FALSE, stringsAsFactors = FALSE)
F_Dicer_KO$sample = "F_Dicer_KO"
F_Dicer_KO = F_Dicer_KO[, 4:5]

F_WT = read.table("mouse_against_zebrafish_STAR_bed_files/F_WT.merged.bed",
                        sep = "\t", header = FALSE, stringsAsFactors = FALSE)
F_WT$sample = "F_WT"
F_WT = F_WT[, 4:5]


F_Adult_WT = read.table("mouse_against_zebrafish_STAR_bed_files/F_Adult_WT.merged.bed",
                        sep = "\t", header = FALSE, stringsAsFactors = FALSE)
F_Adult_WT$sample = "F_Adult_WT"
F_Adult_WT = F_Adult_WT[, 4:5]



# R samples
R_Dicer_KO = read.table("mouse_against_zebrafish_STAR_bed_files/R_Dicer_KO.merged.bed",
                        sep = "\t", header = FALSE, stringsAsFactors = FALSE)
R_Dicer_KO$sample = "R_Dicer_KO"
R_Dicer_KO = R_Dicer_KO[, 4:5]

R_WT = read.table("mouse_against_zebrafish_STAR_bed_files/R_WT.merged.bed",
                  sep = "\t", header = FALSE, stringsAsFactors = FALSE)
R_WT$sample = "R_WT"
R_WT = R_WT[, 4:5]


R_Adult_WT = read.table("mouse_against_zebrafish_STAR_bed_files/R_Adult_WT.merged.bed",
                        sep = "\t", header = FALSE, stringsAsFactors = FALSE)
R_Adult_WT$sample = "R_Adult_WT"
R_Adult_WT = R_Adult_WT[, 4:5]

# Boxplots

cov_counts = rbind(R_WT, R_Dicer_KO, R_Adult_WT, F_WT, F_Dicer_KO, F_Adult_WT)
colnames(cov_counts)[1] = "Coverage_Depth"

cov_counts$sample = factor(cov_counts$sample, levels = c("R_WT", "R_Dicer_KO", "R_Adult_WT", "F_WT", "F_Dicer_KO", "F_Adult_WT"))

p = ggplot(cov_counts, aes(sample, Coverage_Depth, fill = sample)) +
  geom_boxplot() +
  ggtitle("Mouse alignments against Zebrafish", subtitle = "Werner data") + xlab("Sample") + ylab("Depth of Coverage")

p

