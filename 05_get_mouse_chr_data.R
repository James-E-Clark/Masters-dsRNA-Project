# Chromosome barplot for mouse data

library(ggplot2)

server_dir = "http://bsu-srv.ncl.ac.uk/james_clark/primary_alignments/mouse_STAR_bed_files/"
file_suffix = ".primary.reads.per.chr.bed"
column_names = c("chr", "start", "end", "reads")

factor_levels = paste0("chr ", c(1:19, "X", "Y", "M")) # use this vector to order the columns

# F samples
F2095 = read.table(url(paste0(server_dir, "F2095", file_suffix)), stringsAsFactors = FALSE, col.names = column_names)
F2104 = read.table(url(paste0(server_dir, "F2104", file_suffix)), stringsAsFactors = FALSE, col.names = column_names)
F2123 = read.table(url(paste0(server_dir, "F2123", file_suffix)), stringsAsFactors = FALSE, col.names = column_names)
F2128 = read.table(url(paste0(server_dir, "F2128", file_suffix)), stringsAsFactors = FALSE, col.names = column_names)
F2464 = read.table(url(paste0(server_dir, "F2464", file_suffix)), stringsAsFactors = FALSE, col.names = column_names)
F2501 = read.table(url(paste0(server_dir, "F2501", file_suffix)), stringsAsFactors = FALSE, col.names = column_names)
FA1 = read.table(url(paste0(server_dir, "FA1", file_suffix)), stringsAsFactors = FALSE, col.names = column_names)
FA2 = read.table(url(paste0(server_dir, "FA2", file_suffix)), stringsAsFactors = FALSE, col.names = column_names)

# Get table of read counts
F_table = cbind(F2095[, c(1,4)], F2104[, 4], F2123[, 4], F2128[, 4], F2464[, 4], F2501[, 4], FA1[, 4], FA2[, 4])
# multicov counts both reads in the pair, so divide by 2 
F_table[, 2:9] = round(F_table[, 2:9] / 2)
colnames(F_table) = c("chr", "F2095", "F2104", "F2123", "F2128", "F2464", "F2501", "FA1", "FA2")
# Check input proportions against mapping rate
colSums(F_table[, 2:9]) / c(7879473, 7341494, 15490849, 3469568, 10760143, 11694428, 9004878, 7649472)
write.table(F_table, "chromosome_plots/mouse_F_primary_alignments.tsv", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

F_table$chr = factor(gsub("chr", "chr ", F_table$chr), levels = factor_levels)
F_melted = melt(F_table, variable.name = "chr")
colnames(F_melted) = c("chr", "sample", "reads")

g = ggplot(F_melted, aes(chr, reads, fill = sample)) +
  geom_bar(stat = "identity", width = 0.75, position = position_dodge()) +
  ggtitle("Reads per Chromosome (primary alignments)", subtitle = "Werner data") + xlab("Chromosome") + ylab("Total reads") +
  theme(axis.text.x = element_text(angle=65, vjust = 0.5))
g

# Get proportions
F_prop_table = F_table
F_prop_table[, 2:9] = sweep(F_prop_table[, 2:9], 2, colSums(F_prop_table[, 2:9]), '/')

F_prop_table_melted = melt(F_prop_table, variable.name = "chr")
colnames(F_prop_table_melted) = c("chr", "sample", "reads") 

g = ggplot(F_prop_table_melted, aes(chr, reads, fill = sample)) +
  geom_bar(stat = "identity", width = 0.75, position = position_dodge()) +
  ggtitle("Proportion of reads per Chromosome (primary alignments)", subtitle = "Werner data") + xlab("Chromosome") + ylab("Proportion of reads") +
  theme(axis.text.x = element_text(angle=65, vjust = 0.5))
g


# R samples
R2095 = read.table(url(paste0(server_dir, "R2095", file_suffix)), stringsAsFactors = FALSE, col.names = column_names)
R2104 = read.table(url(paste0(server_dir, "R2104", file_suffix)), stringsAsFactors = FALSE, col.names = column_names)
R2123 = read.table(url(paste0(server_dir, "R2123", file_suffix)), stringsAsFactors = FALSE, col.names = column_names)
R2128 = read.table(url(paste0(server_dir, "R2128", file_suffix)), stringsAsFactors = FALSE, col.names = column_names)
R2464 = read.table(url(paste0(server_dir, "R2464", file_suffix)), stringsAsFactors = FALSE, col.names = column_names)
R2501 = read.table(url(paste0(server_dir, "R2501", file_suffix)), stringsAsFactors = FALSE, col.names = column_names)
RA1 = read.table(url(paste0(server_dir, "RA1", file_suffix)), stringsAsFactors = FALSE, col.names = column_names)
RA2 = read.table(url(paste0(server_dir, "RA2", file_suffix)), stringsAsFactors = FALSE, col.names = column_names)

# Get table of read counts
R_table = cbind(R2095[, c(1,4)], R2104[, 4], R2123[, 4], R2128[, 4], R2464[, 4], R2501[, 4], RA1[, 4], RA2[, 4])
# multicov counts both reads in the pair, so divide by 2 
R_table[, 2:9] = round(R_table[, 2:9] / 2)
colnames(R_table) = c("chr", "R2095", "R2104", "R2123", "R2128", "R2464", "R2501", "RA1", "RA2")
# Check input proportions against mapping rate
colSums(R_table[, 2:9]) / c(2339940, 1597443, 2018961, 1999403, 1341157, 2122164, 2168593, 4850284)
write.table(R_table, "chromosome_plots/mouse_R_primary_alignments.tsv", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

R_table$chr = factor(gsub("chr", "chr ", R_table$chr), levels = factor_levels)
R_melted = melt(R_table, variable.name = "chr")
colnames(R_melted) = c("chr", "sample", "reads")

g = ggplot(R_melted, aes(chr, reads, fill = sample)) +
  geom_bar(stat = "identity", width = 0.75, position = position_dodge()) +
  ggtitle("Reads per Chromosome (primary alignments)", subtitle = "Werner data") + xlab("Chromosome") + ylab("Total reads") +
  theme(axis.text.x = element_text(angle=65, vjust = 0.5))
g

# Get proportions
R_prop_table = R_table
R_prop_table[, 2:9] = sweep(R_prop_table[, 2:9], 2, colSums(R_prop_table[, 2:9]), '/')

R_prop_table_melted = melt(R_prop_table, variable.name = "chr")
colnames(R_prop_table_melted) = c("chr", "sample", "reads") 

g = ggplot(R_prop_table_melted, aes(chr, reads, fill = sample)) +
  geom_bar(stat = "identity", width = 0.75, position = position_dodge()) +
  ggtitle("Proportion of reads per Chromosome (primary alignments)", subtitle = "Werner data") + xlab("Chromosome") + ylab("Proportion of reads") +
  theme(axis.text.x = element_text(angle=65, vjust = 0.5))
g

