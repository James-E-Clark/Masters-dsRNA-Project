# ChIPpeakAnno on coverage bedgraphs

setwd("~/Desktop/Gao/")

library(ChIPpeakAnno)
library(org.Mm.eg.db)
library(gtools)
library(EnsDb.Mmusculus.v79)
library(TxDb.Mmusculus.UCSC.mm10.ensGene)

# Function to return a dataframe (bedgraph format) with regions where 
# coverage is >= mult * overall mean coverage  
get_high_coverage_regions = function(bedgraph, mult){
  
  # Get all regions where coverage is greater than mult times the mean   
  peaks = bedgraph[bedgraph$coverage >= mult * mean(bedgraph$coverage),]
  
  # Sort by Chromosome and Start position
  peaks = peaks[order(peaks$space, peaks$start), ]
  
  # Identify contiguous regions
  peaks$contig = numeric(length = nrow(peaks))
  for(i in 2:nrow(peaks)){
    if(paste(peaks$space[i-1], peaks$end[i-1], sep = "_") != paste(peaks$space[i], peaks$start[i], sep = "_")){
      peaks$contig[i] = peaks$contig[i-1] + 1
    } else {
      peaks$contig[i] = peaks$contig[i-1]
    }
  }
  
  # Create data frame to hold results
  region = data.frame(space = character(0), start = numeric(0), end = numeric(0), coverage = numeric(0))
  
  # Get start, end and mean coverage
  for(i in unique(peaks$contig)){
    contig = peaks[peaks$contig == i, ]
    contig$end[1] = contig$end[nrow(contig)]
    contig$coverage[1] = mean(contig$coverage)
    region = rbind(region, contig[1, 1:4])
  }
  
  return(region[mixedorder(region$space), ])
}

#'---------------------
#'Get coverage bedgraph
#'---------------------

bedgraph_merged = read.table("~/Desktop/Gao/INPUT.merged.bed")
colnames(bedgraph_merged) = c("space", "start", "end", "coverage")
cov_bed_merged = get_high_coverage_regions(bedgraph_merged, 3)
write.table(cov_bed_merged, file = "~/Desktop/Gao/Bedgraphs/cov_bed_merged.bed",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Get coverage table with ensembl URLs; make these into hyperlinks in Excel 
location = paste0(gsub("chr", "", cov_bed_merged$space), ":", cov_bed_merged$start, "-", cov_bed_merged$end)
ensembl_url = paste0("https://www.ensembl.org/Mus_musculus/Location/View?r=", location, ";db=core")
(cov_table_merged = cbind(cov_bed_merged, ensembl_url))

#'---------------------------
#'Annotation via ChIPpeakAnno
#'---------------------------

# Prepare annotation data
annoData = toGRanges(EnsDb.Mmusculus.v79, feature = "gene")

# Get GRanges objects for coverage bedgraph file
gr_bed_merged = toGRanges(cov_bed_merged, format = "BED")

gr_bed.anno1 = annotatePeakInBatch(gr_bed_merged, 
                                   AnnotationData=annoData, 
                                   output="nearestBiDirectionalPromoters",
                                   bindingRegion=c(-2000, +2000))


# Summarize the distribution of peaks over different types of features
aCR1 = assignChromosomeRegion(gr_bed_merged,
                              nucleotideLevel = FALSE,
                              precedence = c("Promoters", "immediateDownstream", "fiveUTRs", "threeUTRs", "Exons", "Introns"),
                              TxDb = TxDb.Mmusculus.UCSC.mm10.ensGene)

aCR1
# Plot doesn't fit window, so adjust margins 
par(mar = c(10,4.1,4.1,2.1))
barplot(aCR1$percentage, las = 3)

# Reset margins
par(mar = c(5.1,4.1,4.1,2.1))

# Write to an excel sheet for ease of use
library(WriteXLS)
WriteXLS(cov_table_merged, ExcelFileName = "Gao_INPUT_cov_table.xls", SheetNames = NULL, perl = "perl",

           verbose = FALSE, Encoding = c("UTF-8", "latin1", "cp1252"),

           row.names = FALSE, col.names = TRUE,

           AdjWidth = FALSE, AutoFilter = FALSE, BoldHeaderRow = FALSE,

           na = "",

           FreezeRow = 0, FreezeCol = 0,

           envir = parent.frame())


