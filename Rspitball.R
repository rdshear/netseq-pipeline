# Rspitball.R
suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicRanges)
  library(GenomicAlignments)
  library(rtracklayer)
})

setwd("~/Projects/netseq-pipeline/barcode_tests/umi_tools_from_RX/")

infile <- "dd_200.bam"

raw <- readGAlignments(infile, use.names = TRUE, param = ScanBamParam(tag = "RX"))
bg_pos <- import("~/Projects/netseq-pipeline/barcode_tests/umi_tools_from_RX/dd_200.pos.bedgraph.bgz")
bg_neg <- import("~/Projects/netseq-pipeline/barcode_tests/umi_tools_from_RX/dd_200.neg.bedgraph.bgz")

