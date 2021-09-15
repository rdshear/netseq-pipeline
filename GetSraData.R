library(tidyverse)
library(GEOquery)
# library(GenomicRanges)
# library(GenomicFeatures)
# library(rtracklayer)
# library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
library(SRAdb)

gseID <- "GSE159603"

Sys.setenv("VROOM_CONNECTION_SIZE" = 2^20)
gse <- getGEO(gseID, GSEMatrix = TRUE)[[1]]
p <- phenoData(gse)

sra_dbname <- file.path(system.file('extdata', package='SRAdb'), 'SRAmetadb_demo.sqlite')	
sra_con <- dbConnect(dbDriver("SQLite"), sra_dbname)
