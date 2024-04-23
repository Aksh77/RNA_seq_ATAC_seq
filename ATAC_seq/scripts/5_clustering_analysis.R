library(DiffBind)
library(tidyverse)

# get count matrix
dbObj <- dba(sampleSheet="ATAC_seq/files/sample_sheets/all_samples.csv")
dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE)

# clustering analysis
clustering = "ATAC_seq/data/output_data/clustering_analysis/clustering.png"
png(clustering,width=5,height=5,units="in",res=800, pointsize=8)
plot(dbObj)
dev.off()


