# Get per-chromosome read counts for primary alignments for HeLa, Hilz and Werner data

library(Rsubread)
library(reshape2)
library(ggplot2)

setwd("~/customers/Werner/James_Clark/")

# HeLa data
hela_files = list.files("STAR_HeLa_alignments", pattern = ".primary.bam$", full.names = TRUE)

# Annotation
saf_file = read.table("genome/chr_map_human.tsv", sep = "\t", stringsAsFactors = FALSE, header = TRUE)

# Run featureCounts
map_list = featureCounts(files = hela_files,
                         annot.ext = saf_file,
                         primaryOnly = TRUE,
                         strandSpecific = 0,
                         minMQS = 0,
                         isPairedEnd = TRUE)

# keep a copy of map_list
#save(map_list, file = "HeLa_map_list.RData")
load("HeLa_map_list.RData")

# Extract count table
count_table = map_list$counts
head(count_table)
# Tidy up column names
colnames(count_table) = gsub(".primary.bam", "", colnames(count_table))
count_table = as.data.frame(count_table)

sample_counts = colSums(count_table)
sample_counts

# Compare sample counts with STAR input read numbers
sample_counts / c(17020993, 18492218, 22357752, 18663449, 18076434, 13411906)

# Check rownames are unique
length(unique(rownames(count_table))) == nrow(count_table)

# Get mappings from Ensembl IDs to chromosome
chr_map = unique(saf_file[, 1:2])

# Merge chromosome map with count table
count_table$GeneID = rownames(count_table)

df = merge(count_table, chr_map, by = "GeneID")

# Check that read counts agree with original count table
colSums(df[, 2:7]) == sample_counts

# Get per-chromosome read counts
chr_counts = as.data.frame(colSums(df[df$Chr == "chr1", 2:7]))
colnames(chr_counts) = "chr1"
for(i in paste0("chr", c(2:22, "X", "Y", "M"))){
  chr_counts = cbind(chr_counts, colSums(df[df$Chr == i, 2:7]))
  colnames(chr_counts)[ncol(chr_counts)] = i
}

chr_counts = t(chr_counts)
chr_counts

# Check total counts against STAR input read numbers
colSums(chr_counts) == sample_counts

write.table(chr_counts, "HeLa_STAR_bed_files/HeLa.primary.alignments.per.chr.tsv",
            quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)

# Get barplot of total read counts per chromosome
total_barplot = melt(chr_counts)
colnames(total_barplot) = c("chr", "sample", "reads")

g = ggplot(total_barplot, aes(chr, reads, fill = sample)) +
  geom_bar(stat = "identity", width = 0.75, position = position_dodge()) +
  ggtitle("Reads per chromosome (primary alignments)", subtitle = "HeLa data") + xlab("Chromosome") + ylab("Total reads") +
  theme(axis.text.x = element_text(angle=65, vjust = 0.5))
g

# Get barplot for proportion of reads
prop_barplot = chr_counts[, 1:6] / colSums(chr_counts[, 1:6])
prop_barplot = melt(prop_barplot)
colnames(prop_barplot) = c("chr", "sample", "reads")

g = ggplot(prop_barplot, aes(chr, reads, fill = sample)) +
  geom_bar(stat = "identity", width = 0.75, position = position_dodge()) +
  ggtitle("Proportion of reads per chromosome (primary alignments)", subtitle = "HeLa data") + xlab("Chromosome") + ylab("Proportion of reads") +
  theme(axis.text.x = element_text(angle=65, vjust = 0.5))
g



