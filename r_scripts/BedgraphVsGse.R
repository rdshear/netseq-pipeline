# BedgraphVsGse.R
suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicRanges)
  library(GenomicAlignments)
  library(rtracklayer)
  library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
})

# TODO Parameterize and Cloudify input files
# TODO Debug & un-hard code outliers
new_pos <- import("/n/groups/churchman/rds19/data/S005/wt-1.pos.bedgraph.gz", genome = "sacCer3")
new_neg <- import("/n/groups/churchman/rds19/data/S005/wt-1.neg.bedgraph.gz", genome = "sacCer3")
old_pos <- import("/n/groups/churchman/GSE159603/wt-1.pos.bedgraph.gz", genome = "sacCer3")
old_neg <- import("/n/groups/churchman/GSE159603/wt-1.neg.bedgraph.gz",  genome = "sacCer3")

sacCer3Ranges <- SeqinfoForBSGenome("sacCer3")
tiles <- unlist(tileGenome(sacCer3Ranges, tilewidth = 1000))

binIt <- function(u) {
  result <- binnedAverage(tiles, replace_na(mcolAsRleList(u, "score"), 0), "meanScore", na.rm=TRUE)
  result$normalized <- result$meanScore / mean(result$meanScore)
  result
}

binned_pos_new <- binIt(new_pos)
binned_pos_old <- binIt(old_pos)

result <- cor(binned_pos_new$meanScore, binned_pos_old$meanScore)
print(result)

epsilon <- 1E-9
result <- cor(log(binned_pos_new$meanScore + epsilon), log(binned_pos_old$meanScore + epsilon))
print(result)


v <- data.frame(new = binned_pos_new$meanScore, old = binned_pos_old$meanScore)
plot(v)

filter <- GRanges("chrIV:437850-437860")
old_filtered <- subsetByOverlaps(old_pos, filter)
print(old_filtered)
new_filtered <- subsetByOverlaps(new_pos, filter)
print(new_filtered)

filtered_aligned <- readGAlignments("~/Projects/netseq-pipeline/test_results/wt-1_20210829/outputs/wt-1.aligned.bam", 
                                    param = ScanBamParam(which = filter))

old_wt2_pos <- import("/n/groups/churchman/GSE159603/GSM4835592_wt-2.pos.bedgraph.gz", genome = "sacCer3")
subsetByOverlaps(old_wt2_pos, filter)

bigidx <- (old_wt2_pos$score > 4000)
old_wt2_pos[bigidx]
print(old_wt2_pos[bigidx]$score)
