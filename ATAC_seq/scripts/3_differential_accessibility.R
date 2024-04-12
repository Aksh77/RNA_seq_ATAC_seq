# Load necessary libraries
library(edgeR)
library(EDASeq)
library(GenomicAlignments)
library(GenomicFeatures)
library(wesanderson)
library(Hmisc)
library(dplyr)

ff = FaFile("data/ref/mm10.fa")

# Read count data and define experimental groups
count_file = 'data/output_data/featureCounts/hsc_cfue/hsc_cfue_peaks.counts'
cnt_table = read.table(count_file, sep="\t", header=TRUE, blank.lines.skip=TRUE)
rownames(cnt_table)=cnt_table$Geneid
colnames(cnt_table)=c("Geneid","Chr","Start","End","Strand","Length","HSC1","HSC2","CFUE1","CFUE2")
groups = factor(c(rep("HSC",2),rep("CFUE",2)))
reads.peak = cnt_table[,c(7:10)]

# Get GC content of peak regions for GC-aware normalisation
gr = GRanges(seqnames=cnt_table$Chr, ranges=IRanges(cnt_table$Start, cnt_table$End), strand="*", mcols=data.frame(peakID=cnt_table$Geneid))
peakSeqs = getSeq(x=ff, gr)
gcContentPeaks = letterFrequency(peakSeqs, "GC",as.prob=TRUE)[,1]
gcGroups = Hmisc::cut2(gcContentPeaks, g=20)
mcols(gr)$gc = gcContentPeaks

# Get offsets to correct for library size and GC content
reads.peak=as.matrix(reads.peak)
dataOffset = withinLaneNormalization(reads.peak,y=gcContentPeaks,num.bins=20,which="full",offset=TRUE)
dataOffset = betweenLaneNormalization(reads.peak,which="full",offset=TRUE)

# Calculate differential accessibility using edgeR
design = model.matrix(~groups)
d = DGEList(counts=reads.peak, group=groups)
keep = filterByExpr(d)
d=d[keep,,keep.lib.sizes=FALSE]
d$offset = -dataOffset[keep,]
d.eda = estimateGLMCommonDisp(d, design = design)
fit = glmFit(d.eda, design = design)
lrt.EDASeq = glmLRT(fit, coef = 2)
DA_res=as.data.frame(topTags(lrt.EDASeq, nrow(lrt.EDASeq$table)))
DA_res$Geneid = rownames(DA_res)
DA.res.coords = left_join(DA_res,cnt_table[1:4],by="Geneid")

# Filter for significant peaks with FDR < 0.05
DA.res.coords = DA.res.coords[DA_res$FDR < 0.05,]

# Save results
outdir = "data/output_data/differential_accessibility"
if (file.exists(outdir) == FALSE) {
  dir.create(outdir,recursive=TRUE)
}
setwd(outdir)
result_file = "hsc_vs_cfue.tsv"
write.table(DA.res.coords, result_file, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, fileEncoding = "")

# make volcano plot
pdf("hsc_vs_cfue_volcano.pdf", width=5, height=5)
plot(DA_res$logFC, -log10(DA_res$FDR), pch=20, col=ifelse(DA_res$FDR < 0.05, "red", "black"),
     xlab="log2 fold change", ylab="-log10 FDR", main="HSC vs CFUE differential accessibility")
dev.off()

# make an MA plot
pdf("hsc_vs_cfue_MA.pdf", width=5, height=5)
plot(DA_res$logCPM, DA_res$logFC, pch=20, col=ifelse(DA_res$FDR < 0.05, "red", "black"),
     xlab="Average log2 counts", ylab="log2 fold change", main="HSC vs CFUE differential accessibility")
dev.off()
