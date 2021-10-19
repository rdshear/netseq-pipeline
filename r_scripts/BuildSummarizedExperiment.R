# BuildRaggedExperiment.R
# Create RaggedExperiment object for NETseq screen data set
# Add derived statistics that are functions of the scores and possibly the changepoints
library(GenomicRanges)
library(rtracklayer)
library(SummarizedExperiment)
library(fitdistrplus)
library(tidyverse)
library(tidygenomics)
library(plyranges)

set.seed(20190413)


dfile.path <- "/n/groups/churchman/rds19/data/S005/"
outfile.path <- "/Users/robertshear/Documents/n/groups/churchman/rds19/data/S005/HresSE_screen.rds"
feature.filename <- str_c(dfile.path, "genelist.gff")
cpa_algorithm <- "CEZINB"
max_genes <- 8

# apply function b to each element of matrix a, returning a matrix of shape a
mat.sapply <- function(a, b) {
  structure(matrix(sapply(a, b), nrow = nrow(a), ncol = ncol(a)),
            dimnames = dimnames(a))
}

bedgraph_to_granges <- function(pos, neg)
{
  scores <- mapply(function(strand_sym, infilename) {
    x <- import(infilename, genome = "sacCer3")
    strand(x) <- strand_sym
    x
  }, list('+', '-'), list(pos, neg), SIMPLIFY = FALSE)
  scores <- do.call(c, scores)
  x <- gaps(scores)
  x <- x[as.character(strand(x)) != "*"]
  x$score <- 0
  scores <- sort(c(scores, x))
  w <- findOverlaps(feature_ranges, scores)
  ov <- split(subjectHits(w), queryHits(w))
  m <- lapply(seq_along(ov), function(i) {
    u <- restrict(scores[ov[[i]]], start = start(feature_ranges[i]), 
                  end = end(feature_ranges[i]))
    GPos(feature_ranges[i], score = rep(u$score, width(u)))
  })
  matrix(m, ncol = 1)
}

gff3_to_granges <- function(fname) {
  x <- import(fname)
  # Keep only the mcols needed
  mcols(x) <- mcols(x)[,c("tx_name", "seq_index", "m", "v")]
  x$seq_index <- as.integer(x$seq_index)
  x$m <- as.numeric(x$m)
  x$v <- as.numeric(x$v)
  x <- split(x, x$tx_name)
  matrix(lapply(names(feature_ranges), function(u) x[u][[1]]), ncol = 1)
}

features <- import(feature.filename)
# TODO: TEST ONLY
if (max_genes > 0) {
  features <- features[1:max_genes]
}
names(features) <- features$ID
feature_ranges <- GRanges(as_tibble(features)[,c("seqnames", "start", "end", "strand")])
names(feature_ranges) <- names(features)

sample_table <- tibble(sample_id = c("wt-1", "wt-2")) %>%
    mutate(variant = map_chr(sample_id,
                         function(u) str_split(u, "-")[[1]][1]),
    bedgraph_neg = str_c(dfile.path, sample_id, ".neg.bedgraph.gz"),
    bedgraph_pos = str_c(dfile.path, sample_id, ".pos.bedgraph.gz"),
    changepoints = str_c(dfile.path, sample_id, ".cp.", cpa_algorithm, ".gff3")) %>%
    column_to_rownames(var = "sample_id")

e <- do.call(cbind, lapply(seq(nrow(sample_table)), function(i) {
  u <- sample_table[i,]
  x <- gff3_to_granges(u$changepoints)
  y <- bedgraph_to_granges(u$bedgraph_pos, u$bedgraph_neg)
  r <- SummarizedExperiment(rowData = features, colData = u)
  assays(r, withDimnames = FALSE) <- SimpleList(segments = x, scores = y)
  r
  }))

assay(e, "k") <- mat.sapply(assay(e, "segments"), length)

e1 <- e[8,1]
segs <- assay(e1, "segments")[[1]]
counts <- assay(e1, "scores")[[1]]
s <- counts$score
# TODO WRAP THESE to MATs

  x <- try(fitdistrplus::fitdist(s, "nbinom"), silent = TRUE)
  if (inherits(x, "try-error")) {
    alpha <- NA;
    alpha.se <- NA;
    mu.se <- NA
  } else {
    alpha <- x$estimate["size"]
    alpha.se <- x$sd["size"]
    mu.se <- x$sd["mu"]
  }

  if (length(segs) < 2) {
    cliff.magnitude <- NA
  } else {
    if (as.character(strand(segs[1])) == "-") {
      x <- rev(segs)
    } else {
      x <- segs
    }
    cliff.magnitude <- mcols(x[2])$m / mcols(x[1])$m
  }
#   
#   list(cp = cp, mu = mean(s), var = var(s), 
#        k = length(cp), 
#        alpha = alpha, 
#        alpha.se = alpha.se, 
#        mu.se = mu.se,
#        cliff.magnitude = cliff.magnitude
#        )
#   }, 
#   cp = segments(e), s = occupancyRle(e))
# 
# for (assay.name in rownames(new.assays)) {
#   assay(e, assay.name) <- reroll(new.assays[assay.name, ], e)
# }
# 
# # fixup superflouous list for scalar assays
# for (u in names(assays(e))) {
#   if (class(assay(e, u)[[1]]) %in% c("numeric", "integer")) {
#     x <- assay(e, u)
#     assay(e, u) <- reroll(unlist(unroll(x)), x)
#   }
# }

saveRDS(e, file = outfile.path)
#
