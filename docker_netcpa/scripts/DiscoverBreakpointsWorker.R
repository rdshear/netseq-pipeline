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
  setwd("~/temp")
  print("DEBUG IS ON -- COMMAND LINE PARAMETERS IGNORED")
  commandArgs <- function(trailingOnly) {
    c("~/temp/shard_1.rds",
      "~/temp/cp_shard_1.gff",
      "50", # Maximum gene body length
      "12" # Kmax (maximum number of segments)
      )
  }
}
#  .............TO HERE

args <- commandArgs(trailingOnly = TRUE)

input.filename <- args[1]
output.filename <- args[2]
GeneMaxLength <- as.numeric(args[3]) # Truncate gene to this length (or inf if 0)
Kmax <- as.numeric(args[4])

source_data <- readRDS(args[1])

print(sprintf("Starting at %s.  Shard name = %s. Genes = %d", 
        Sys.time(), input.filename, length(source_data)))

algorithm <- "CEZINB"


options(mc.cores = detectCores())
print(sprintf("Number of cores detected = %d", getOption("mc.cores")))


start.time <- Sys.time()

result <- mclapply(source_data, function(u) {
  gene <- u$gene
  # truncate gene lengths if so desired
  if (GeneMaxLength > 0) {
    gene <- resize(gene, fix = "start", ifelse(width(gene) > GeneMaxLength, GeneMaxLength, width(gene)))
  }
  
  s <- u$scores
  # if there is no overlap between the feature (u) and the bedGraph entries,
  # then mcolAsRleList will fail. Workaround follows
  if (length(s) == 0) {
    s <- Rle(values = 0, lengths = width(u))
  } else {
    s <- mcolAsRleList(s, "score")[[seqnames(gene)]][start(gene):end(gene)]
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
  
  result <- GRanges(seqnames = rep(seqnames(gene)[1], n),
              strand = rep(strand(gene)[1], n),
              ranges = IRanges(start = start(gene) - 1 + tr[, 1], 
                end = start(gene) - 1 + tr[, 2]),
              tx_name = gene$ID,
              seq_index = seq(n),
              type = "seq_index",
              source = "DiscoverBP",
             algorithm = rep(algorithm, n),
             m = round(stats["m", ],4),
             v = round(stats["v", ],4))
  
  
  
    writeLines(sprintf("TRACE %s  Gene %s\n", Sys.time(), gene$ID), stderr())
    writeLines(kableExtra::kable(gc(), format = "pipe"), stderr())
    close(fileConn)
  
    result
  },
  # mclapply paramters
  mc.preschedule = FALSE, mc.silent = FALSE, mc.cores = 1)

result <- unlist(GRangesList(result))

end.time <- Sys.time()
elapsed.time <- difftime(end.time, start.time, units = "secs")
time.per.gene <- elapsed.time / length(source_data)
total_width <- sum(sapply(source_data, function(u) min(width(u$gene), GeneMaxLength)))
time.per.nt <- elapsed.time / total_width

  
  # # TODO: Carry BIC and logLikelihood
  # # # TODO: get segment statistics
  # # # TODO: reverse negative strand 
  # # IRanges(start = c(1, seg + 1), end = c(seg, length(s)))
  # # # TODO: report multi-map removal areas
  # c(tx_name = u$tx_name, mu = mu, v = v, 
  #   mzl = head(sort(runLength(s), decreasing = TRUE)))

export.gff3(result, con = output.filename)

cat(sprintf("Elapsed time: %.0f sec,  genes: %.0f,   bases: %.0f \n  %0.2f sec/gene,   %0.1f msec/base \n Completed at %s\n",
        elapsed.time, length(source_data), total_width, time.per.gene, time.per.nt * 1000, Sys.time()))

print(sessionInfo())
