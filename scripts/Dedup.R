library(tidyverse)
library(GenomicRanges)
library(GenomicAlignments)
library(rtracklayer)

wd <- "/Users/robertshear/temp/outputs/"
infile <- "wt-1.aligned.bam"

range <- "chrVI:1-1000000"
sample <- "wt-1"

wd <- "/Users/robertshear/temp/outputs/"
infile <- "xwt-1.aligned.bam"

range <- NULL
sample <- "xwt-1"


timestamp(suffix = " Start")

setwd(wd)

readGAlignments(infile, use.names = FALSE, param = ScanBamParam(tag = "RX", which= GRanges(range))) %>%
  print %>%
  keep(qwidth(.) - width(.) > 0) %>%
  granges(use.mcols=TRUE) %>%
  resize(fix="start", ignore.strand=FALSE, width=1) %>%
  as_tibble %>%
  group_by(seqnames, start, strand) %>%
  summarise(score = n_distinct(RX)) %>%
  mutate(end = start) %>% 
  GRanges -> result
strand(result) <- if_else(as.vector(strand(result)) == '+', '-','+')
timestamp(suffix = " Read Complete")

for (strand in c("+", "-")) {
  outfile = paste0(sample, ".", ifelse(strand == "+", "pos","neg"), ".bedgraph")
  export.bedGraph(result[as.character(strand(result)) == strand], outfile, index=TRUE)
  timestamp(suffix = paste("", outfile, " Export Complete"))
}
