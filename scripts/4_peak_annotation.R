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

# Annotate peaks with closest genomic features
bed.annot = annotatePeak(peaks.gr, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Mm.eg.db")
annot_peaks = as.data.frame(bed.annot)

# Save to file
outdir = "data/output_data/peak_annotation"
if (file.exists(outdir) == FALSE) {
  dir.create(outdir,recursive=TRUE)
}
setwd(outdir)
result_file = "hsc_vs_cfue.tsv"
write.table(annot_peaks, result_file,
            append = FALSE,
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            fileEncoding = "")

# Visualise the annotation summary
pdf("hsc_vs_cfue.pdf")
upsetplot(bed.annot, vennpie=TRUE)
dev.off()
