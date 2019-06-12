# ggplot barplot example

library(ggplot2)

# Example: Hilz data, barplot wih proportion of uniquely mapped reads for day1, day3 and day7

# Could use functions or loops here (especially at the outset), but doing it long-hand for clarity

server_dir = "http://bsu-srv.ncl.ac.uk/james_clark/Hilz_Bowtie2_bed_files/"
file_suffix = ".unique.reads.per.chr.bed"
column_names = c("chr", "start", "end", "reads")

factor_levels = paste0("chr ", c(1:19, "X", "Y", "M")) # use this vector to order the columns

# Read in the data from bsu-srv
day1_n1 = read.table(url(paste0(server_dir, "day1_whole_testis_n1", file_suffix)), stringsAsFactors = FALSE, col.names = column_names)
day1_n2 = read.table(url(paste0(server_dir, "day1_whole_testis_n2", file_suffix)), stringsAsFactors = FALSE, col.names = column_names)

day3_n1 = read.table(url(paste0(server_dir, "day3_whole_testis_n1", file_suffix)), stringsAsFactors = FALSE, col.names = column_names)
day3_n2 = read.table(url(paste0(server_dir, "day3_whole_testis_n2", file_suffix)), stringsAsFactors = FALSE, col.names = column_names)

day7_n1 = read.table(url(paste0(server_dir, "day7_whole_testis_n1", file_suffix)), stringsAsFactors = FALSE, col.names = column_names)
day7_n2 = read.table(url(paste0(server_dir, "day7_whole_testis_n2", file_suffix)), stringsAsFactors = FALSE, col.names = column_names)

# Merge and get total reads for each day
day1_df = merge(day1_n1[, c(1,4)], day1_n2[, c(1,4)], by = "chr")
day1_df$total = day1_df$reads.x + day1_df$reads.y # getting total read count
day1_df$proportion = day1_df$total / sum(day1_df$total) # getting read counts as a proportion of total
day1_df$day = "day_1" # adding a column which ggplot will use to separate bars in barplot

day3_df = merge(day3_n1[, c(1,4)], day3_n2[, c(1,4)], by = "chr")
day3_df$total = day3_df$reads.x + day3_df$reads.y
day3_df$proportion = day3_df$total / sum(day3_df$total)
day3_df$day = "day_3"

day7_df = merge(day7_n1[, c(1,4)], day7_n2[, c(1,4)], by = "chr")
day7_df$total = day7_df$reads.x + day7_df$reads.y
day7_df$proportion = day7_df$total / sum(day7_df$total)
day7_df$day = "day_7"

# ggplot needs a 'long' data frame with 'chr', 'total' and 'day' columns
# rbind stacks data frames vertically
total_barplot = rbind(day1_df[, c(1,4,6)], day3_df[, c(1,4,6)], day7_df[, c(1,4,6)])
# Set the order of chromosomes (default is a-z, 0-9)
total_barplot$chr = gsub("chr", "chr ", total_barplot$chr) # adding a space for readability on plot
total_barplot$chr = factor(total_barplot$chr, levels = factor_levels)

g = ggplot(total_barplot, aes(chr, total, fill = day)) +
  geom_bar(stat = "identity", width = 0.75, position = position_dodge()) +
  ggtitle("Total reads per Chromosome", subtitle = "Hilz data") + xlab("Chromosome") + ylab("Total reads") +
  theme(axis.text.x = element_text(angle=65, vjust = 0.5))
g

# To use proportion of reads, use 'chr', 'proportion' and 'day' columns
prop_barplot = rbind(day1_df[, c(1,5,6)], day3_df[, c(1,5,6)], day7_df[, c(1,5,6)])
prop_barplot$chr = gsub("chr", "chr ", prop_barplot$chr)
prop_barplot$chr = factor(prop_barplot$chr, levels = factor_levels)

g2 = ggplot(prop_barplot, aes(chr, proportion, fill = day)) +
  geom_bar(stat = "identity", width = 0.75, position = position_dodge()) +
  ggtitle("Proportion of uniquely mapped reads per Chromosome", subtitle = "Hilz data") + xlab("Chromosome") + ylab("Proportion of reads") +
  theme(axis.text.x = element_text(angle=65, vjust = 0.5))
g2


