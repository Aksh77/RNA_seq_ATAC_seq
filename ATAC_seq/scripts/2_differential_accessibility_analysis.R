# Load necessary libraries
library(DiffBind)
library(tidyverse)
library(ChIPseeker)
library(GenomicAlignments)
library(GenomicFeatures)
library(clusterProfiler)
library(biomaRt)
library(ReactomePA)
library(ggupset)
library(ggimage)
library(EnhancedVolcano)
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# read contrast from command line arguments
# Usage: Rscript 4_peak_annotation_and_functional_enrichment.R hsc_cfue
args = commandArgs(trailingOnly=TRUE)
contrast = args[1]
contrast = "hsc_cfue"

# ****************** FIND DIFFERENTIAL ACCESSIBILE REGIONS ********************

# make directory for saving differential accessibility results
diff_acc_outdir = paste("ATAC_seq/data/output_data/differential_accessibility/",contrast, sep="")
if (file.exists(diff_acc_outdir) == FALSE) {
  dir.create(diff_acc_outdir,recursive=TRUE)
}

# read in peaksets
sample_sheet = paste("ATAC_seq/files/sample_sheets/",contrast,".csv", sep="")
dbObj <- dba(sampleSheet=sample_sheet)

# get count information for each of the peaks
dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE)

# Heatmap of samples
heatmap = paste(diff_acc_outdir,"/",contrast,"_heatmap.png", sep="")
png(heatmap,width=5,height=5,units="in",res=800, pointsize=8)
plot(dbObj)
dev.off()

# establish contrasts
dbObj <- dba.contrast(dbObj, categories=DBA_FACTOR, minMembers = 2)

# differential enrichment analysis
dbObj <- dba.analyze(dbObj, method=DBA_ALL_METHODS)
dba.show(dbObj, bContrasts=T)

# Venn diagram for overlap of significant regions between methods
venn_plot = paste(diff_acc_outdir,"/",contrast,"_venn.png", sep="")
png(venn_plot,width=10,height=10,units="in",res=1200, pointsize=20)
dba.plotVenn(dbObj,contrast=1,method=DBA_ALL_METHODS)
dev.off()

# Volcano plot
volcano_plot = paste(diff_acc_outdir,"/",contrast,"_volcano.png", sep="")
png(volcano_plot,width=5,height=5,units="in",res=800, pointsize=8)
dba.plotVolcano(dbObj)
dev.off()


# MA plot
ma_plot = paste(diff_acc_outdir,"/",contrast,"_MA.png", sep="")
png(ma_plot,width=5,height=5,units="in",res=800, pointsize=12)
dba.plotMA(dbObj, contrast=1)
dev.off()

# save list of significant differential regions
res <- dba.report(dbObj, contrast=1, th=1)
res <- res[res$FDR < 0.05,]
diff_acc_results = paste("ATAC_seq/data/output_data/differential_accessibility/",contrast,"/",contrast,"_diff_acc.csv", sep="")
write.csv(res, file=diff_acc_results)

# # ****************************** PEAK ANNOTATION ******************************

# make directory for saving peak annotation results
outdir = paste("ATAC_seq/data/output_data/peak_annotation/",contrast, sep="")
if (file.exists(outdir) == FALSE) {
  dir.create(outdir,recursive=TRUE)
}

# Annotate peaks in in differentially accessible regions with closest genomic featuresdiff_acc_file =
diff_acc_file = paste("ATAC_seq/data/output_data/differential_accessibility/",contrast,"/",contrast,"_diff_acc.csv", sep="")
res = read.table(diff_acc_file, sep=",", header=TRUE)

# peak annotation
peaks.gr = GRanges(seqnames=res$seqnames, ranges=IRanges(res$start, res$end), strand=res$strand)
bed.annot = annotatePeak(peaks.gr, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Mm.eg.db")
annot_peaks = as.data.frame(bed.annot)

# plot pie chart of genomic features
feature_plot = paste(outdir,"/",contrast,"_genomic_features.png", sep="")
png(feature_plot,width=5,height=5,units="in",res=800, pointsize=8)
plotAnnoPie(bed.annot)
dev.off()

# merge differential accessibility and peak annotation using seqnames, start and end
annot_peaks$PeakID = paste(annot_peaks$seqnames, annot_peaks$start, annot_peaks$end, sep="_")
res = as.data.frame(res)
res$PeakID = paste(res$seqnames, res$start, res$end, sep="_")
annot_peaks = merge(annot_peaks, res, by="PeakID")

# plot enhanced volcano plot
labels = annot_peaks$SYMBOL
enhanced_volcano_plot = paste(diff_acc_outdir,"/",contrast,"_enhanced_volcano.png", sep="")
png(enhanced_volcano_plot,width=7,height=7,units="in",res=800, pointsize=8)
EnhancedVolcano(annot_peaks, x = 'Fold', y = 'FDR', lab=labels)
dev.off()

# Save to file
result_file = paste(outdir,"/",contrast,"_annotated_peaks.tsv", sep="")
write.table(annot_peaks, result_file,
            append = FALSE,
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            fileEncoding = "")

# Visualise the annotation summary
summary_plot = paste(outdir,"/",contrast,"_annotated_peaks.pdf", sep="")
pdf(summary_plot)
upsetplot(bed.annot, vennpie=TRUE)
dev.off()

# Plot distribution of peaks relative to TSS
tss_dist_plot = paste(outdir,"/",contrast,"_dist_to_tss.pdf", sep="")
pdf(tss_dist_plot)
plotDistToTSS(bed.annot, title="Distribution of ATAC-seq peaks loci\nrelative to TSS")
dev.off()

# *********************** GENE SET ENRICHMENT ANALYSIS ************************

# make directory for saving GO enrichment results
func_outdir = paste("ATAC_seq/data/output_data/GO_enrichment/",contrast, sep="")
if (file.exists(func_outdir) == FALSE) {
  dir.create(func_outdir,recursive=TRUE)
}

# get top 500 most significantly differentially accessible upregulated and downregulated peaks 
annot_peaks$Fold = as.numeric(annot_peaks$Fold)
annot_up = annot_peaks[annot_peaks$Fold > 0,]
upregulated_peaks = annot_up[order(annot_up$FDR),][1:500,]
annot_down = annot_peaks[annot_peaks$Fold < 0,]
downregulated_peaks = annot_down[order(annot_down$FDR),][1:500,]

# Find enriched pathways
pathway.reac <- enrichPathway(as.data.frame(top_annot_peaks)$geneId, organism="mouse")
result_file = paste(func_outdir,"/",contrast,"_enriched_pathways.tsv", sep="")
write.table(pathway.reac, result_file,
            append = FALSE,
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            fileEncoding = "")

# Find enriched GO terms for Molecular Function
upregulated_mf <- enrichGO(as.data.frame(upregulated_peaks)$geneId, org.Mm.eg.db, ont = "MF")
result_file = paste(func_outdir,"/",contrast,"_enriched_GO_MF_upregulated.tsv", sep="")
write.table(upregulated_mf, result_file,
            append = FALSE,
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            fileEncoding = "")
downregulated_mf = enrichGO(as.data.frame(downregulated_peaks)$geneId, org.Mm.eg.db, ont = "MF")
result_file = paste(func_outdir,"/",contrast,"_enriched_GO_MF_downregulated.tsv", sep="")
write.table(downregulated_mf, result_file,
            append = FALSE,
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            fileEncoding = "")

# Find enriched GO terms for Biological Process
upregulated_bp <- enrichGO(as.data.frame(upregulated_peaks)$geneId, org.Mm.eg.db, ont = "BP")
result_file = paste(func_outdir,"/",contrast,"_enriched_GO_BP_upregulated.tsv", sep="")
write.table(upregulated_bp, result_file,
            append = FALSE,
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            fileEncoding = "")
downregulated_bp = enrichGO(as.data.frame(downregulated_peaks)$geneId, org.Mm.eg.db, ont = "BP")
result_file = paste(func_outdir,"/",contrast,"_enriched_GO_BP_downregulated.tsv", sep="")
write.table(downregulated_bp, result_file,
            append = FALSE,
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            fileEncoding = "")

