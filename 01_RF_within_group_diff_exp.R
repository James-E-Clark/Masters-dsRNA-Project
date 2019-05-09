#-------------------------------------------------------------------------------------------#
#  Author      : John Casement                                                              |
#  Language    : R Statistical Programming Language                                         |
#  Data Owner  : Andi Werner                                                                |
#  System      : Sulaco                                                                     |
#  Tool        : R - DESeq2                                                                 |
#  Description : Run DESeq2 on Salmon output                                                |
#-------------------------------------------------------------------------------------------#

setwd("/data/customers/Natasya/Salmon_analysis/")

library(DESeq2)
library(readr)
library(tximport)
library(biomaRt)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(gtools)
library(annotate)


##'Read in Transcript to Gene Map
##'-----------------------------------------------------------------------------------------#
tx2gene <- read_tsv("mm10/mm10_gene_map.tsv", col_names = c("txid", "geneid"))
##'-----------------------------------------------------------------------------------------#


##'Read in Salmon Abundances, and set up DESeq2 Dataset - EXCLUDING KO samples
##'-----------------------------------------------------------------------------------------#

# N.B. fastq at john@beatles:~/customers/Natasya/fastq

files_in = mixedsort(list.files(".", recursive = TRUE, pattern = "*.sf"))
files_in

# Remove outlier-ish samples R2464 and F2104
files_in = files_in[grep("R2464", files_in, invert = TRUE)]
files_in = files_in[grep("F2104", files_in, invert = TRUE)]
files_in

sample_names = gsub("/quant.sf", "", gsub("Salmon/", "", files_in))
sample_names

# set up pheno data

# Dicer KO samples are 2123, 2104 and 2464
# WT samples are 2501, 2128 and 2095 
# Adult WT samples are A1 and A2

pheno = data.frame(sample_name = sample_names,
                   group = NA,
                   type = substr(sample_names, 1, 1))

for(i in c("2123","2104","2464")) pheno$group[grep(i, pheno$sample_name)] = "Dicer"
for(i in c("2501","2128","2095")) pheno$group[grep(i, pheno$sample_name)] = "WT"
for(i in c("A1", "A2")) pheno$group[grep(i, pheno$sample_name)] = "Adult_WT"

# Split by F and R groups
pheno_F = pheno[pheno$type == "F",]
pheno_R = pheno[pheno$type == "R",]

files_F = files_in[substr(files_in, 8, 8) == "F"]
files_R = files_in[substr(files_in, 8, 8) == "R"]

#  read in Salmon abundances
txi_F = tximport(files   = files_F, 
                 type    = "salmon", 
                 tx2gene = tx2gene)

txi_R = tximport(files   = files_R, 
                 type    = "salmon", 
                 tx2gene = tx2gene)

# create DESeqDataSets for between outcomes comparison 
deseq_dataset_F = DESeqDataSetFromTximport(txi     = txi_F,
                                           colData = pheno_F, 
                                           design  = ~ group)

deseq_dataset_R = DESeqDataSetFromTximport(txi     = txi_R,
                                           colData = pheno_R, 
                                           design  = ~ group)

# check colData correspondence
colData(deseq_dataset_F)
colData(deseq_dataset_R)

# run DESeq
deseq_dataset_F = DESeq(deseq_dataset_F)
colnames(deseq_dataset_F) = as.vector(colData(deseq_dataset_F)$sample_name)

deseq_dataset_R = DESeq(deseq_dataset_R)
colnames(deseq_dataset_R) = as.vector(colData(deseq_dataset_R)$sample_name)

# get raw counts and normalised counts
count_table_F = counts(deseq_dataset_F, normalized = FALSE)
head(count_table_F, 10)
write.table(count_table_F, file = "F_samples/count_table_F.tsv", sep ="\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
count_table_norm_F = counts(deseq_dataset_F, normalized = TRUE)
head(count_table_norm_F, 10)
write.table(count_table_norm_F, file = "F_samples/count_table_normalized_F.tsv", sep ="\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

count_table_R = counts(deseq_dataset_R, normalized = FALSE)
head(count_table_R, 10)
write.table(count_table_R, file = "R_samples/count_table_R.tsv", sep ="\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
count_table_norm_R = counts(deseq_dataset_R, normalized = TRUE)
head(count_table_norm_R, 10)
write.table(count_table_norm_R, file = "R_samples/count_table_normalized_F.tsv", sep ="\t", row.names = TRUE, col.names = TRUE, quote = FALSE)


##'-----------------------------------------------------------------------------------------#
##'Sample distances and PCA: F samples
##'-------------------------------------------------------------------------------------------#

library("RColorBrewer")
library("gplots")
library("pheatmap")

colors = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

# Sample distances
rlog_data_F = rlog(deseq_dataset_F, blind=TRUE)

sampleDists_F = dist(t(assay(rlog_data_F)))
sampleDistMatrix_F = as.matrix(sampleDists_F)
colnames(sampleDistMatrix_F) = paste(pheno_F$group)
rownames(sampleDistMatrix_F) = paste(pheno_F$group)

pheatmap(sampleDistMatrix_F,
         clustering_distance_rows=sampleDists_F,
         clustering_distance_cols=sampleDists_F,
         col=colors)

# PCA
library("ggplot2")
library("ggrepel")
plot_data = plotPCA(rlog_data_F, intgroup="group", returnData=TRUE)
percentVar = round(100 * attr(plot_data, "percentVar"))
ggplot(plot_data, aes(PC1, PC2, color=group, label = name)) +
#  geom_label_repel(aes(PC1, PC2, color=group, label = name)) +
  geom_text(hjust = "left", nudge_x = 0.2) +
  geom_point(size=2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

##'-------------------------------------------------------------------------------------------#


##'-----------------------------------------------------------------------------------------#
##'Sample distances and PCA: R samples
##'-------------------------------------------------------------------------------------------#

library("RColorBrewer")
library("gplots")
library("pheatmap")

colors = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

# Sample distances
rlog_data_R = rlog(deseq_dataset_R, blind=TRUE)

sampleDists_R = dist(t(assay(rlog_data_R)))
sampleDistMatrix_R = as.matrix(sampleDists_R)
colnames(sampleDistMatrix_R) = paste(pheno_R$group)
rownames(sampleDistMatrix_R) = paste(pheno_R$group)

pheatmap(sampleDistMatrix_R,
         clustering_distance_rows=sampleDists_R,
         clustering_distance_cols=sampleDists_R,
         col=colors)

# PCA
plot_data = plotPCA(rlog_data_R, intgroup="group", returnData=TRUE)
percentVar = round(100 * attr(plot_data, "percentVar"))
ggplot(plot_data, aes(PC1, PC2, color=group, label = name)) +
  #  geom_label_repel(aes(PC1, PC2, color=group, label = name)) +
  geom_text(hjust = "left", nudge_x = 0.2) +
  geom_point(size=2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

##'-------------------------------------------------------------------------------------------#


##' Get annotation
##'-------------------------------------------------------------------------------------------# 

all(rownames(count_table_norm_F) == rownames(count_table_R)) # TRUE

library(biomaRt)

ensembl = useMart("ensembl")
ensembl = useDataset("mmusculus_gene_ensembl",
                     mart = ensembl)

annot = getBM(attributes = c('ensembl_gene_id',  'mgi_symbol',
                             'description',      'chromosome_name',
                             'start_position',   'end_position',
                             'strand',           'gene_biotype'),
              filters    = 'ensembl_gene_id',
              values     = rownames(count_table_F),
              mart       = ensembl)

#write.table(annot, file = "mmusculus_gene_annot.tsv", sep ="\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
annot = read.table("mmusculus_gene_annot.tsv", sep = "\t", row.names = 1, header = TRUE, stringsAsFactors = FALSE, quote = "\"")

##'-------------------------------------------------------------------------------------------#


##'-----------------------------------------------------------------------------------------#
##'Get results
##'-------------------------------------------------------------------------------------------#

resultsNames(deseq_dataset_F)
resultsNames(deseq_dataset_R)

run_diff_exp = function(dds, group_1, group_2, pval_cut, fc_cut){

  res = results(dds, contrast = c("group", group_1, group_2))
  
  # remove NAs
  res = res[complete.cases(res),]
  # include ensembl_id column (for merge with annot)
  res$ensembl_gene_id = rownames(res)
  # add annotation
  res_anno = merge(as.data.frame(res), annot, by = "ensembl_gene_id")
  # sort by magnitude of fold change
  res_anno = res_anno[rev(order(abs(res_anno$log2FoldChange))),]
  # get significant genes by log fold change and adjusted p-value and write to file
  res_sig = res_anno[(res_anno$padj < pval_cut) & (abs(res_anno$log2FoldChange) > log2(fc_cut)), ]
  
  return(list(res_sig = res_sig, res_anno = res_anno))

}

# F samples
df = run_diff_exp(dds = deseq_dataset_F,  group_1 = "Adult_WT", group_2 = "WT", pval_cut = 0.05, fc_cut = 2)
dim(df$res_sig) # 59 14
write.table(df$res_sig, file = "F_samples/Adult_WT_vs_WT.sig.genes.tsv" , sep = "\t", row.names = FALSE, quote = FALSE)
write.table(df$res_anno, file = "F_samples/Adult_WT_vs_WT.all.genes.tsv" , sep = "\t", row.names = FALSE, quote = FALSE)

df = run_diff_exp(dds = deseq_dataset_F,  group_1 = "Dicer", group_2 = "WT", pval_cut = 0.05, fc_cut = 2)
dim(df$res_sig) # nothing passes cut-offs
#write.table(df$res_sig, file = "F_samples/Dicer_vs_WT.sig.genes.tsv" , sep = "\t", row.names = FALSE, quote = FALSE)
write.table(df$res_anno, file = "F_samples/Dicer_vs_WT.all.genes.tsv" , sep = "\t", row.names = FALSE, quote = FALSE)

# R samples
df = run_diff_exp(dds = deseq_dataset_R,  group_1 = "Adult_WT", group_2 = "WT", pval_cut = 0.05, fc_cut = 2)
dim(df$res_sig) # 206 14
write.table(df$res_sig, file = "R_samples/Adult_WT_vs_WT.sig.genes.tsv" , sep = "\t", row.names = FALSE, quote = FALSE)
write.table(df$res_anno, file = "R_samples/Adult_WT_vs_WT.all.genes.tsv" , sep = "\t", row.names = FALSE, quote = FALSE)

df = run_diff_exp(dds = deseq_dataset_R,  group_1 = "Dicer", group_2 = "WT", pval_cut = 0.05, fc_cut = 2)
dim(df$res_sig) # 41 14
write.table(df$res_sig, file = "R_samples/Dicer_vs_WT.sig.genes.tsv" , sep = "\t", row.names = FALSE, quote = FALSE)
write.table(df$res_anno, file = "R_samples/Dicer_vs_WT.all.genes.tsv" , sep = "\t", row.names = FALSE, quote = FALSE)


##'-------------------------------------------------------------------------------------------#

library(ChIPpeakAnno)

