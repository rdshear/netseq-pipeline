# DiscoverBreakpointsWorker.R
# bedgraph files --> gff3 files
# appropriate for scatter/gather model
# 
# Usage:
# Rscript --vanilla scripts/DiscoverBreakpointsWorker.R \
#   {shard_n.rds}    # changepoint input
#   {changepoints.gff}    # changepoint output
#   {K.max}   # maximum number of changepoints
#

# To run with embedded parameters, set DEBUG.TEST <- TRUE

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(breakpoint)
  library(tidyverse)
  library(plyranges)
  library(parallel)
  library(foreach)
  library(doMC)
})
set.seed(20210915)

# DEBUG ONLY FROM HERE.....
DEBUG.TEST <- TRUE
if (interactive() && exists("DEBUG.TEST")) {
  setwd("~/temp")
  print("DEBUG IS ON -- COMMAND LINE PARAMETERS IGNORED")
  commandArgs <- function(trailingOnly) {
    c("~/temp/shard_0.rds",
      "~/temp/cp_shard_0.gff",
      "12" # Kmax (maximum number of segments)
      )
  }
}
#  .............TO HERE

args <- commandArgs(trailingOnly = TRUE)

input.filename <- args[1]
output.filename <- args[2]
Kmax <- as.numeric(args[3])

print(sprintf("Starting at %s.  Shard name = %s.", 
        Sys.time(), input.filename))

algorithm <- "CEZINB"

get_change_points <- function(s) {
  u <- data.frame(score = s)
  # NOTE: This is *always* + direction. We assume no bias in this agorithm by strand
  seg <- try(CE.ZINB(data = u, Nmax = Kmax, parallel = FALSE), silent = FALSE)
  if (inherits(seg, c("character", "try-error"))) {
    tau <- numeric(0)
  } else {
    tau <- seg$BP.Loc
  }
  result <- tibble(seq_index = seq(length(tau) + 1),
         s.start = c(1, tau), s.end = c(tau - 1, nrow(u))) %>%
    mutate(m = map2_dbl(s.start, s.end,
                        function(a, b, c) mean(c[a:b]), u$score),
                      #TODO {var} or {v}
                      #TODO Round m and v
                      var = map2_dbl(s.start, s.end,
                          function(a, b, c) var(c[a:b]), u$score))
  result
}

mc.cores <- max(detectCores(), 2) - 1
registerDoMC(cores = mc.cores)
start.time <- Sys.time()
print(sprintf("Starting at %s, Number of cores to use = %d",
             as.character(start.time), mc.cores))

readRDS(input.filename) %>% 
  mutate(scores = pmap(list(start, end, data), function(s, e, d) {
      locs <- map2(d$start - s + 1, d$end - d$start + 1, function(u, v) u + seq(0, v-1))
      x <- rep(0, e - s + 1)
      # generate the score. loop sometimes beats obsucrity
      for (i in seq_along(locs)) {
        x <- modify_at(x, locs[[i]], function(a,b) a + b, d$score[i])
      }
      x
      })) %>%
  mutate(mu = map_dbl(scores, mean), v = map_dbl(scores, var)) %>%
  mutate(segments = foreach(s = .$scores) %dopar% get_change_points(s)) %>%
  unnest(segments) -> result 


gr.out <- GRanges(seqnames = result$seqnames,
                  strand = result$strand,
                  ranges = IRanges(start = result$s.start, end = result$s.end),
                  tx_name = result$ID,
                  seq_index = result$seq_index,
                  type = "seq_index",
                  algorithm = algorithm,
                  m = round(result$m,4),
                  v = round(result$var,4))

end.time <- Sys.time()
elapsed.time <- difftime(end.time, start.time, units = "secs")
gene.count <- length(unique(result$ID))
time.per.gene <- elapsed.time / gene.count
total_width <- sum(result$end - result$start)
time.per.nt <- elapsed.time / total_width

  
  # # TODO: Carry BIC and logLikelihood
  # # # TODO: get segment statistics
  # # # TODO: reverse negative strand 
  # # IRanges(start = c(1, seg + 1), end = c(seg, length(s)))
  # # # TODO: report multi-map removal areas
  # c(tx_name = u$tx_name, mu = mu, v = v, 
  #   mzl = head(sort(runLength(s), decreasing = TRUE)))

export.gff3(gr.out, con = output.filename)

cat(sprintf("Elapsed time: %.0f sec,  genes: %.0f,   bases: %.0f \n  %0.2f sec/gene,   %0.1f msec/base \n Completed at %s\n",
        elapsed.time, gene.count, total_width, time.per.gene, time.per.nt * 1000, Sys.time()))

print(sessionInfo())
