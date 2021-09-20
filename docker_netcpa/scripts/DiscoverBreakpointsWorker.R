# DiscoverBreakpointsWorker.R
# bedgraph files --> gff3 files
# appropriate for scatter/gather model
# 
# Usage:
# Rscript --vanilla scripts/DiscoverBreakpointsWorker.R \
#   {regions.of.ineterest.gff} #regions of interest
#   {occupancy.bedgraph.pos.gz} 
#   {occupancy.bedgraph.neg.gz}
#   {changepoints.gff}    # changepoint output
#   {max.gene.length}     # regions of interest longer than this will only 
#                         # be searched to this lengh
#   {K.max}   # maximum number of changepoints
#

# To run with embedded parameters, set DEBUG.TEST <- TRUE

suppressPackageStartupMessages({
  library(parallel)
  library(GenomicRanges)
  library(rtracklayer)
  library(breakpoint)
})
set.seed(20210915)

# DEBUG ONLY FROM HERE.....
DEBUG.TEST <- TRUE
if (interactive() && exists("DEBUG.TEST")) {
  print("DEBUG IS ON -- COMMAND LINE PARAMETERS IGNORED")
  commandArgs <- function(trailingOnly) {
    c("~/temp/shard_1.gff",
      "/n/groups/churchman/rds19/data/S005/wt-1.pos.bedgraph.gz",
      "/n/groups/churchman/rds19/data/S005/wt-1.neg.bedgraph.gz",
      "~/temp/cp_shard_1.gff",
      "300", # Maximum gene body length
      "12" # Kmax (maximum number of segments)
      )
  }
}
#  .............TO HERE

args <- commandArgs(trailingOnly = TRUE)

subject_genes.filename <- args[1]
infile.pos <- args[2]
infile.neg <- args[3]
output.filename <- args[4]
GeneMaxLength <- as.numeric(args[5]) # Truncate gene to this length (or inf if 0)
Kmax <- as.numeric(args[6])

g <- import(subject_genes.filename)

sprintf("Starting at %s.  Shard name = %s. Genes = %d", 
        Sys.time(), subject_genes.filename, length(g))

algorithm <- "CEZINB"


options(mc.cores = detectCores())
sprintf("Number of cores detected = %d", getOption("mc.cores"))

names(g) <- g$ID

# truncate gene lengths if so desired
if (GeneMaxLength > 0) {
  g <- resize(g, fix = "start", ifelse(width(g) > GeneMaxLength, GeneMaxLength, width(g)))
}

start.time <- Sys.time()

result <- mclapply(as(g, "GRangesList"), function(u) {
  is.plus <- as.logical(as.character(strand(u)) == "+")
  s <- import(ifelse(is.plus, infile.pos, infile.neg), which = u, genome = "sacCer3")
  # if there is no overlap between the feature (u) and the bedGraph entries,
  # then mcolAsRleList will fail. Workaround follows
  if (length(s) == 0) {
    s <- Rle(values = 0, lengths = width(u))
  } else {
    s <- mcolAsRleList(s, "score")[[seqnames(u)]][start(u):end(u)]
  }
  s <- replace(s, is.na(s), 0)

  mu <- mean(s)
  v <- var(s)

  # NOTE: This is *always* + direction. We assume no bias in this agorithm by strand
  seg <- try(CE.ZINB(data = data.frame(s), Nmax = Kmax, parallel = FALSE), silent = TRUE)
  if (inherits(seg, c("character", "try-error"))) {
    tau <- numeric(0)
  } else {
    tau <- seg$BP.Loc
  }
  n <- length(tau) + 1
  
  tr <- cbind(c(1, tau), c(tau - 1, length(s)))

  
  stats <- apply(tr, 1, function(w) {
    v <- s[w[1]:w[2]]
    c(m = mean(v), v = var(v))
  })
  
  result <- GRanges(seqnames = rep(seqnames(u)[1], n),
              strand = rep(strand(u)[1], n),
              ranges = IRanges(start = start(u) - 1 + tr[, 1], 
                end = start(u) - 1 + tr[, 2]),
              seq_index = seq(n),
              type = "seq_index",
              source = "DiscoverBP",
             algorithm = rep(algorithm, n),
             m = round(stats["m", ],4),
             v = round(stats["v", ],4))
  result
})

# HACK. Can't set the name in the GRanges constructor

for (i in seq_along(result)) {
  result[[i]]$tx_name = names(result)[i]
}

result <- unlist(GRangesList(result))

end.time <- Sys.time()
elapsed.time <- difftime(end.time, start.time, units = "secs")
time.per.gene <- elapsed.time / length(g)
time.per.nt <- elapsed.time / sum(width(g))

#  mcols(result)$tx_name <- as.character(u$tx_name)
  #result$type <- Rle(values = "seq_index", lengths = n)
  
  # # TODO: Carry BIC and logLikelihood
  # # # TODO: get segment statistics
  # # # TODO: reverse negative strand 
  # # IRanges(start = c(1, seg + 1), end = c(seg, length(s)))
  # # # TODO: report multi-map removal areas
  # c(tx_name = u$tx_name, mu = mu, v = v, 
  #   mzl = head(sort(runLength(s), decreasing = TRUE)))

export.gff3(result, con = output.filename)

cat(sprintf("Elapsed time: %.0f sec,  genes: %.0f ,   bases: %.0f \n  %0.2f sec/gene,   %0.1f msec/base \n Completed at %s\n",
        elapsed.time, length(g), sum(width(g)), time.per.gene, time.per.nt * 1000, Sys.time()))

print(sessionInfo())
