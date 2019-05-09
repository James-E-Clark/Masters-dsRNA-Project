# ChIPpeakAnno on mouse mm10 coverage bedgraphs

setwd("~/customers/Werner/James_Clark/")

library(ChIPpeakAnno)
library(org.Mm.eg.db)
library(EnsDb.Mmusculus.v79)
library(gtools)

# Function to return a dataframe (bedgraph format) with regions where coverage is >= mult * overall mean coverage  
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

sample_name = "R2095"

sample_URL = paste0("http://bsu-srv.ncl.ac.uk/james_clark/bed_files/", sample_name, ".bedgraph")

sample_bedgraph = read.table(url(sample_URL))
colnames(sample_bedgraph) = c("space", "start", "end", "coverage")

# Examples:

# Regions where coverage is at least 8 times the overall mean 
cov_bed = get_high_coverage_regions(sample_bedgraph, 8)
# Regions where coverage is at least 3 times the overall mean 
cov_bed = get_high_coverage_regions(sample_bedgraph, 2)

# Get coverage table with ensembl URLs; you could make these into hyperlinks in Excel 
location = paste0(gsub("chr", "", cov_bed$space), ":", cov_bed$start, "-", cov_bed$end)
ensembl_url = paste0("https://www.ensembl.org/Mus_musculus/Location/View?r=", location, ";db=core")
(cov_table = cbind(cov_bed, ensembl_url))


#'---------------------------
#'Annotation via ChIPpeakAnno
#'---------------------------

# Prepare annotation data
annoData = toGRanges(EnsDb.Mmusculus.v79, feature = "gene")

# Get GRanges objects for coverage bedgraph file
gr_bed = toGRanges(cov_bed, format = "BED")

gr_bed.anno = annotatePeakInBatch(gr_bed, 
                                  AnnotationData=annoData, 
                                  output="nearestBiDirectionalPromoters",
                                  bindingRegion=c(-2000, +2000))
gr_bed.anno

# Summarize the distribution of peaks over different types of features
aCR = assignChromosomeRegion(gr_bed,
                             nucleotideLevel = FALSE,
                             precedence = c("Promoters", "immediateDownstream", "fiveUTRs", "threeUTRs", "Exons", "Introns"),
                             TxDb = TxDb.Mmusculus.UCSC.mm10.ensGene)

aCR

# Plot doesn't fit window, so adjust margins 
par(mar = c(10,4.1,4.1,2.1))
barplot(aCR$percentage, las = 3)
# Reset margins
par(mar = c(5.1,4.1,4.1,2.1))
