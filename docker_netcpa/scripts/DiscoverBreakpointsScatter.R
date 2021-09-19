# DiscoverBreakpointsScatter.R
# This step distributes the reference gene ranges to
# an array of files, each of which will have ranges for the 
# worker task to execute in gff3 format
# 
# Usage:
# Rscript --vanilla scripts/DiscoverBreakpointsScatter.R \
#   {reference.gene.list.filename} \
#   {n.genes} \
#   {n.shards} \ 
#   {output.file.directory}
#
# Example:
# Rscript --vanilla scripts/DiscoverBreakpointsScatter.R \
    # /n/groups/churchman/rds19/data/S005/genelist.gff \
    # 12 \
    # 2 \
    # /n/groups/churchman/rds19/data/S005/ 

# To run with embedded parameters, set DEBUG.TEST <- TRUE

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
})
set.seed(20190416)

# DEBUG ONLY FROM HERE.....
DEBUG.TEST <- TRUE
if (interactive() && exists("DEBUG.TEST")) {
  print("DEBUG IS ON -- COMMAND LINE PARAMETERS IGNORED")
  commandArgs <- function(trailingOnly) {
    c("/n/groups/churchman/rds19/data/S005/genelist.gff",
      "5", # n.genes
      "2", # n.shards
      "~/temp/scatter_list/")
  }
}
#  .............TO HERE

args <- commandArgs(trailingOnly = TRUE)

subject_genes.filename <- args[1]
maxGenes <- as.numeric(args[2]) # if > 0 sample this number of genes
n.shards <- as.numeric(args[3]) # shards
output.directory <- args[4]


sprintf("Starting at %s.  shards = %s. Max genes = %d", 
        Sys.time(), n.shards, maxGenes)

g <- import(subject_genes.filename)

# subset the genes if so desired
if (maxGenes > 0 & maxGenes < length(g)) {
  g <- g[sort(sample(length(g), maxGenes))]
}

# create the shards
gs <- split(g, rep_len(seq(1, n.shards), length(g)))

# TODO Remove unnneded mcols
dir.create(output.directory, recursive = TRUE)
for (i in seq_along(gs)) {
  fn <- file.path(output.directory, paste0("shard_",i , ".gff"))
  export.gff3(gs[[i]], fn)
}

print(sessionInfo())
