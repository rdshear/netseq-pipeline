# Pool Samples
library(rtracklayer)
library(glue)
# TODO: INVESTIGATE Lis Lab Package!
library(BRGenomics)

dir <- "/n/groups/churchman/rds19/data/S005"
pat <- glue("^wt.*[.]dedup[.]bam$")



ps_list <- lapply( list.files(path = dir, pattern = pat, full.names = TRUE), function(bfile) {
  bfile <- "/n/groups/churchman/rds19/data/S005/wt-1.dedup.bam"
  ps <- import_bam(bfile, 
                   mapq = 20, 
                   revcomp = TRUE,
                   trim.to = "3p",
                   paired_end = FALSE)
  ps <- ps[as.character(seqnames(ps)) != "chrM"]
  print(bfile)
  print(length(ps))
  ps
})

ps_pool <- mergeGRangesData(ps_list)

strands <- c(pos = "+", neg = "-")
for (i in seq_along(strands)) {
  outfile_name <- glue("{dir}/wt-pool.{names(strands[i])}.bedgraph.gz")
  x <- ps_pool[as.character(strand(ps_pool)) == strands[i]]
  print(outfile_name)
  print(length(x))
  export.bedGraph(x, outfile_name)
}
