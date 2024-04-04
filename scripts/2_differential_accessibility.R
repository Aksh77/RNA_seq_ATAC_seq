library(optparse)
library(DESeq2)
library(vsn)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(lattice)
library(BiocParallel)

# Read counts table and format it
count.table <- read.delim(file='data/output_data/hsc_cfue_peaks.counts',header=TRUE,skip=1)
colnames(count.table) <- gsub(".bam","",colnames(count.table))
colnames(count.table) <- as.character(lapply(colnames(count.table), function (x) tail(strsplit(x,'.',fixed=TRUE)[[1]],1)))
rownames(count.table) <- count.table$Geneid
interval.table <- count.table[,1:6]
count.table <- count.table[,7:ncol(count.table),drop=FALSE]

# Create output directory for DEseq2
outdir='data/output_data/hsc_cfue'
if (file.exists(outdir) == FALSE) {
  dir.create(outdir,recursive=TRUE)
}
setwd(outdir)

# Run DEseq2
samples.vec <- sort(colnames(count.table))
groups <- sub("_[^_]+$", "", samples.vec)
print(unique(groups))

# outprefix='hsc_cfue'
# DDSFile <- paste(outprefix,".dds.rld.RData",sep="")
# if (file.exists(DDSFile) == FALSE) {
#   counts <- count.table[,samples.vec,drop=FALSE]
#   coldata <- data.frame(row.names=colnames(counts),condition=groups)
#   dds <- DESeqDataSetFromMatrix(countData = round(counts), colData = coldata, design = ~ condition)
#   dds <- DESeq(dds, parallel=TRUE, BPPARAM=MulticoreParam(opt$cores))
#   if (!opt$vst) {
#     rld <- rlog(dds)
#   } else {
#     rld <- vst(dds)
#   }
#   save(dds,rld,file=DDSFile)
}

