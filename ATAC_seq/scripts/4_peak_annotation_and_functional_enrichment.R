# Load necessary libraries
library(ChIPseeker)
library(GenomicAlignments)
library(GenomicFeatures)
library(clusterProfiler)
library(biomaRt)
library(ReactomePA)
library(ggupset)
library(ggimage)
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# read contrast from command line arguments
# Usage: Rscript 4_peak_annotation_and_functional_enrichment.R hsc_cfue
args = commandArgs(trailingOnly=TRUE)
contrast = args[1]

# make directory for saving peak annotation results
outdir = paste("ATAC_seq/data/output_data/peak_annotation/",contrast, sep="")
if (file.exists(outdir) == FALSE) {
  dir.create(outdir,recursive=TRUE)
}

# make directory for saving functional enrichment results
func_outdir = paste("ATAC_seq/data/output_data/functional_enrichment/",contrast, sep="")
if (file.exists(func_outdir) == FALSE) {
  dir.create(func_outdir,recursive=TRUE)
}

# Annotate peaks in in differentially accessible regions with closest genomic featuresdiff_acc_file =
diff_acc_file = paste("ATAC_seq/data/output_data/differential_accessibility/",contrast,"/",contrast,"_differential_accessibility.tsv", sep="")
DA.res.coords = read.table(diff_acc_file, sep="\t", header=TRUE)
peaks.gr = GRanges(seqnames=DA.res.coords$Chr, ranges=IRanges(DA.res.coords$Start, DA.res.coords$End), strand=DA.res.coords$Strand)
bed.annot = annotatePeak(peaks.gr, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Mm.eg.db")
annot_peaks = as.data.frame(bed.annot)

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

# Find enriched pathways
pathway.reac <- enrichPathway(as.data.frame(annot_peaks)$geneId, organism="mouse")
pathway.reac[1:10,c(1:7,9)]

# Find enriched GO terms for Molecular Function
pathway.GO <- enrichGO(as.data.frame(annot_peaks)$geneId, org.Mm.eg.db, ont = "MF")

# Save results
result_file = paste(func_outdir,"/",contrast,"_enriched_GO.tsv", sep="")
write.table(pathway.reac, result_file,
            append = FALSE,
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            fileEncoding = "")
