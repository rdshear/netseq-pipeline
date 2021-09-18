# DiscoverBreakpointsWorker.R
# 
# Usage:
# Rscript --vanilla scripts/DiscoverBreakpointsGather.R \
#   {changepoints.gff}    # combined output
#   {shard1.gff} {shard2.gff} ... {shardn.gff}
#

# To run with embedded parameters, set DEBUG.TEST <- TRUE

suppressPackageStartupMessages({
  library(parallel)
  library(GenomicRanges)
  library(rtracklayer)
  library(breakpoint)
})
set.seed(20210916)

# DEBUG ONLY FROM HERE.....
DEBUG.TEST <- TRUE
if (interactive() && exists("DEBUG.TEST")) {
  print("DEBUG IS ON -- COMMAND LINE PARAMETERS IGNORED")
  commandArgs <- function(trailingOnly) {
    c("~/temp/scatter_list/wt-1.cp.gff",
      "~/temp/scatter_list/shard_1.cp.gff",
      "~/temp/scatter_list/shard_2.cp.gff"
      )
  }
}
#  .............TO HERE

args <- commandArgs(trailingOnly = TRUE)

output.filename <- args[1]
shards <- args[-1]

result <- unlist(GRangesList(lapply(shards, function(u) import(u))))

# TODO: add seqinfo and sort


export(result, con = output.filename)

sprintf("Completed at %s\n", Sys.time())
print(sessionInfo())
