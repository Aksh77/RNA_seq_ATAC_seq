library(DiffBind)
library(tidyverse)

dbObj <- dba(sampleSheet='files/samples.csv')

# get consensus peaks
dbObj.cons <- dba.peakset(dbObj, consensus = DBA_FACTOR, minOverlap = 2)
dbObj.cons <- dba(dbObj.cons, mask = dbObj.cons$masks$Consensus, minOverlap = 1)
dba.plotVenn(dbObj.cons, mask = 1:2)
consensus.peaks <- dba.peakset(dbObj.cons, bRetrieve = TRUE)

dbObj.counts <- dba.count(dbObj, peaks = consensus.peaks)

# visualising replicates
dba.plotPCA(dbObj,  attributes=DBA_FACTOR, label=DBA_ID)
plot(dbObj)

dbObj <- dba.contrast(dbObj, categories=DBA_FACTOR, minMembers = 2)
