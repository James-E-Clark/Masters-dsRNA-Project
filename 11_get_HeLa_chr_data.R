# Chromosome barplot for HeLa data

library(ggplot2)

server_dir = "http://bsu-srv.ncl.ac.uk/james_clark/primary_alignments/HeLa_STAR_bed_files/"
file_suffix = ".primary.reads.per.chr.bed"
column_names = c("chr", "start", "end", "reads")

factor_levels = paste0("chr ", c(1:22, "X", "Y", "M"))

# Read in the data from bsu-srv
untreated = read.table(url(paste0(server_dir, "Untreated", file_suffix)), stringsAsFactors = FALSE, col.names = column_names)
untreated = untreated[, c(1,4)]
untreated$proportion = untreated$reads / sum(untreated$reads)
untreated$sample = "Untreated"

siCntrl = read.table(url(paste0(server_dir, "siCntrl", file_suffix)), stringsAsFactors = FALSE, col.names = column_names)
siCntrl = siCntrl[, c(1,4)]
siCntrl$proportion = siCntrl$reads / sum(siCntrl$reads)
siCntrl$sample = "siCntrl"

siSUV3_rep1 = read.table(url(paste0(server_dir, "siSUV3_rep1", file_suffix)), stringsAsFactors = FALSE, col.names = column_names)
siSUV3_rep1 = siSUV3_rep1[, c(1,4)]
siSUV3_rep1$proportion = siSUV3_rep1$reads / sum(siSUV3_rep1$reads)
siSUV3_rep1$sample = "siSUV3_rep1"

siSUV3_rep2 = read.table(url(paste0(server_dir, "siSUV3_rep2", file_suffix)), stringsAsFactors = FALSE, col.names = column_names)
siSUV3_rep2 = siSUV3_rep2[, c(1,4)]
siSUV3_rep2$proportion = siSUV3_rep2$reads / sum(siSUV3_rep2$reads)
siSUV3_rep2$sample = "siSUV3_rep2"

siPNPase_rep1 = read.table(url(paste0(server_dir, "siPNPase_rep1", file_suffix)), stringsAsFactors = FALSE, col.names = column_names)
siPNPase_rep1 = siPNPase_rep1[, c(1,4)]
siPNPase_rep1$proportion = siPNPase_rep1$reads / sum(siPNPase_rep1$reads)
siPNPase_rep1$sample = "siPNPase_rep1"

siPNPase_rep2 = read.table(url(paste0(server_dir, "siPNPase_rep2", file_suffix)), stringsAsFactors = FALSE, col.names = column_names)
siPNPase_rep2 = siPNPase_rep2[, c(1,4)]
siPNPase_rep2$proportion = siPNPase_rep2$reads / sum(siPNPase_rep2$reads)
siPNPase_rep2$sample = "siPNPase_rep2"

# Get table of read counts
HeLa_table = cbind(untreated[, c(1,2)], siCntrl[, 2], siSUV3_rep1[, 2], siSUV3_rep2[, 2], siPNPase_rep1[, 2], siPNPase_rep2[, 2])
# multicov counts both reads in the pair, so divide by 2 
HeLa_table[, 2:7] = round(HeLa_table[, 2:7] / 2)
colnames(HeLa_table) = c("chr", "untreated", "siCntrl", "siSUV3_rep1", "siSUV3_rep2", "siPNPase_rep1", "siPNPase_rep2")
write.table(HeLa_table, "chromosome_plots/HeLa_primary_alignments.tsv", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

# make long data frame for ggplot
total_barplot = rbind(untreated, siCntrl, siSUV3_rep1, siSUV3_rep2, siPNPase_rep1, siPNPase_rep2)

# Set the order of chromosomes (default is a-z, 0-9)
total_barplot$chr = gsub("chr", "chr ", total_barplot$chr)
total_barplot$chr = factor(total_barplot$chr, levels = factor_levels)

g = ggplot(total_barplot, aes(chr, reads, fill = sample)) +
  geom_bar(stat = "identity", width = 0.75, position = position_dodge()) +
  ggtitle("Reads per Chromosome (primary alignments)", subtitle = "HeLa data") + xlab("Chromosome") + ylab("Total reads") +
  theme(axis.text.x = element_text(angle=65, vjust = 0.5))
g

# To use proportion of reads, use 'chr', 'proportion' and 'day' columns
prop_barplot = rbind(untreated[, c(1,3,4)], siCntrl[, c(1,3,4)], siSUV3_rep1[, c(1,3,4)], siSUV3_rep2[, c(1,3,4)], siPNPase_rep1[, c(1,3,4)], siPNPase_rep2[, c(1,3,4)])
prop_barplot$chr = gsub("chr", "chr ", prop_barplot$chr)
prop_barplot$chr = factor(prop_barplot$chr, levels = factor_levels)

g2 = ggplot(prop_barplot, aes(chr, proportion, fill = sample)) +
  geom_bar(stat = "identity", width = 0.75, position = position_dodge()) +
  ggtitle("Proportion of reads per Chromosome (primary alignments)", subtitle = "HeLa data") + xlab("Chromosome") + ylab("Proportion of reads") +
  theme(axis.text.x = element_text(angle=65, vjust = 0.5))
g2

