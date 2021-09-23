# DiscoverBreakpointsScatter.R
# This step distributes the reference gene ranges to
# an array of files, each of which will have ranges for the 
# worker task to execute in gff3 format
# 
# Usage:
# Rscript --vanilla scripts/DiscoverBreakpointsScatter.R \
#   {reference.gene.list.filename} \
#   {occupancy.bedgraph.pos.gz} 
#   {occupancy.bedgraph.neg.gz}
#   {n.genes} \
#   {n.shards} 
#
# Example:
# Rscript --vanilla scripts/DiscoverBreakpointsScatter.R \
    # /n/groups/churchman/rds19/data/S005/genelist.gff \
    # 12 \
    # 2

# To run with embedded parameters, set DEBUG.TEST <- TRUE

suppressPackageStartupMessages({
  library(rtracklayer)
})
set.seed(20190416)

# DEBUG ONLY FROM HERE.....
DEBUG.TEST <- TRUE
if (interactive() && exists("DEBUG.TEST")) {
  print("DEBUG IS ON -- COMMAND LINE PARAMETERS IGNORED")
  setwd("~/temp/")
  commandArgs <- function(trailingOnly) {
    c("/n/groups/churchman/rds19/data/S005/genelist.gff",
      "/n/groups/churchman/rds19/data/S005/wt-1.pos.bedgraph.gz",
      "/n/groups/churchman/rds19/data/S005/wt-1.neg.bedgraph.gz",
      "5", # n.genes
      "2") # n.shards
  }
}
#  .............TO HERE

args <- commandArgs(trailingOnly = TRUE)

subject_genes.filename <- args[1]
bedgraph.filename.pos <- args[2]
bedgraph.filename.neg <- args[3]
maxGenes <- as.numeric(args[4]) # if > 0 sample this number of genes
n.shards <- as.numeric(args[5]) # shards


sprintf("Starting at %s.  shards = %s. Max genes = %d", 
        Sys.time(), n.shards, maxGenes)

g <- import.gff3(subject_genes.filename, genome = "sacCer3", 
                  feature.type = "gene", colnames = "ID")

x <- mapply(function (filename, strand) {
                  result <- import(filename, genome = "sacCer3")
                  strand(result) <- strand
                  result
                },
              c(bedgraph.filename.pos, bedgraph.filename.neg), c('+', '-'))

scores <- c(x[[1]], x[[2]])

# subset the genes if so desired
if (maxGenes > 0 & maxGenes < length(g)) {
  g <- g[sort(sample(length(g), maxGenes))]
}

v <- findOverlaps(g, scores)
v <- split(subjectHits(v), queryHits(v))
gindex <- as.integer(names(v))
names(gindex) <- names(g[gindex])
w <- mapply(function(gene, idx) {
          list(gene = g[gene], scores = scores[idx])
        },
    gindex, v, SIMPLIFY = FALSE)


# create the shards
gs <- split(w, rep_len(seq(1, n.shards), length(w)))

for (i in seq_along(gs)) {
  fn <- file.path(paste0("shard_",i , ".rds"))
  saveRDS(gs[[i]], fn)
}

print(sessionInfo())
