# Rscript --vanilla /scripts/Dedup.R <output-directory> <input-bam> <sample-name> [<GRanges string>]
library(tidyverse)
library(GenomicRanges)
library(GenomicAlignments)
library(rtracklayer)

rm(list = ls())
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  args <- c("/Users/robertshear/Projects/netseq-pipeline/test/wd/", "deduplicated.bam", "wt-1")
}
wd <- args[1]
infile <- args[2]
sample <- args[3]
range <- args[4]
if (is.na(range)) range <- NULL

timestamp(suffix = " Start")

setwd(wd)

raw <- readGAlignments(infile, use.names = FALSE, param = ScanBamParam(tag = "RX", which= GRanges(range)))

raw %>% keep(njunc(.) & as.character(seqnames(.)) != "chrM") %>%
  granges(use.mcols=TRUE) %>%
  resize(fix="start", ignore.strand=FALSE, width=1) %>%
  as_tibble %>%
  group_by(seqnames, start, strand) %>%
  summarise(score = n()) %>%
  mutate(end = start) %>% 
  GRanges -> result

strand(result) <- if_else(as.vector(strand(result)) == '+', '-','+')
timestamp(suffix = " Read Complete")

for (strand in c("+", "-")) {
  outfile = paste0(sample, ".", ifelse(strand == "+", "pos","neg"), ".bedgraph")
  export.bedGraph(result[as.character(strand(result)) == strand], outfile, index=TRUE)
  timestamp(suffix = paste("", outfile, " Export Complete"))
}

