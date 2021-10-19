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
feature.filename <- str_c(dfile.path, "genelist.gff")
cpa_algorithm <- "CEZINB"

bedgraph_to_granges <- function(pos, neg)
{
  scores <- mapply(function(strand_sym, infilename) {
    x <- import(infilename, genome = "sacCer3")
    strand(x) <- strand_sym
    x
  }, list('+', '-'), list(pos, neg), SIMPLIFY = FALSE)
  sort(unlist(GRangesList(scores)))
}

# TODO Assumes there are no missing changepoint regions of interest
gff3_to_granges <- function(fname) {
  x <- import(fname)
  x <- split(x, x$tx_name)
  as(x[names(feature_ranges)], "SimpleList")
}

# TODO: Store segments as vector of widths
# s <- r[99,1];s
# t <- rowRanges(s);t
# v <- x[[99]];v
# w <- width(v);w
# y <-IRanges(start = start(t) + c(0,cumsum(w))[1:length(w)], width = w);y
# assay(s, "z", withDimnames = FALSE) <- matrix(IntegerList(w), ncol = 1, nrow = 1)
# assay(s, "z")


features <- import(feature.filename)
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

# mutate(...
# gr_post = map2(bedgraph_pos, bedgraph_neg, bedgraph_to_granges),
# gr_cpa = map(changepoints, gff3_to_granges)

u <- sample_table[1,]
x <- gff3_to_granges(u$changepoints)
r <- SummarizedExperiment(rowData = features, colData = u)
assay(r, "changepoint", withDimnames = FALSE) <- matrix(x, ncol = 1)


    # HACK Conflict with tidyverse because GRanges doesnt support '[[' etc


q <- map(sample_table$changepoints, gff3_to_granges)


x <- gff3_to_granges(u$changepoints)
assay(exp, "segments") <- x

z <- transmute(x, sample_id, variant, gr_cpa, gr_post)
# nest the Granges within the gene
z <- RaggedExperiment(colData = sample_table, rowData = features)

# create data.frame of samples
x <- config$samples
u <- x[[1]]
y <- lapply(x, function(u) data.frame(variant = u$variant, group = u$group, stringsAsFactors = FALSE))
z <- Reduce(rbind, y, data.frame())
z$sample <- names(y)
rownames(z) <- z$sample

# TODO: For test only
z <- z[z$group == sample.group, ]

g <- import(feature.filename)
names(g) <- g$ID
# TODO: For test only
# g <- sample(g, 32)
 
e <- HresSE(sample = HresSamples(scoreFileDirectory = occupancy.path,
                                 segmentFileDirectory = changepoint.path,
                                 sampleNames = z$sample, 
                                 scoreFilesPos = paste0(z$sample, ".pos.bedgraph.gz"), 
                                 scoreFilesNeg = paste0(z$sample, ".neg.bedgraph.gz"), 
                                 segmentFiles = paste0(z$sample, ".cp.CEZINB.gff3")),
                    rowRanges = g)

cd <- colData(e)
cd$variant <- z$variant
cd$group <- z$group
colData(e) <- cd

# TODO: THIS IS TEMPORARY UNTIL WE GET THE side-effects figured out
e <- readOccupancy(e)
e <- readSegments(e)

unroll <- function(m) sapply(m, identity)

reroll <- function(v, m)   structure(matrix(v, nrow = nrow(m), ncol = ncol(m)),
                                     dimnames = dimnames(m))

mat.sapply <- function(a, b) {
  structure(matrix(sapply(a, b), nrow = nrow(a), ncol = ncol(a)),
            dimnames = dimnames(a))
}


new.assays <- mapply(function(cp, s) {
  # compute the mean and variance for each segment
  cpx <- shift(cp, 1 - min(start(cp)))
  x <- apply(cbind(start(cpx), end(cpx)), 1, function(u) 
    (function(x) list(mu = mean(x), var = var(x)))(s[u[1]:u[2]]))
  mcols(cp) <- structure(data.frame(matrix(unlist(x), nrow = length(cp), byrow = TRUE)), 
                         names = c("mu", "var"))
  # TODO: change try to try/catch and capture and dispose of the parameter to fitdist's "catch"
  x <- try(fitdistrplus::fitdist(as.numeric(s), "nbinom"), silent = TRUE)
  if (inherits(x, "try-error")) {
    alpha <- NA;
    alpha.se <- NA;
    mu.se <- NA
  } else {
    alpha <- x$estimate["size"]
    alpha.se <- x$sd["size"]
    mu.se <- x$sd["mu"]
  }
  
  # TODO: use hrseSE methods to simplify this
  if (length(cp) < 2) {
    cliff.magnitude <- 0
  } else {
    if (as.character(strand(cp[1])) == "-") {
      x <- rev(cp)
    } else {
      x <- cp
    }
    cliff.magnitude <- mcols(x[2])$mu / mcols(x[1])$mu
  }
  
  list(cp = cp, mu = mean(s), var = var(s), 
       k = length(cp), 
       alpha = alpha, 
       alpha.se = alpha.se, 
       mu.se = mu.se,
       cliff.magnitude = cliff.magnitude
       )
  }, 
  cp = segments(e), s = occupancyRle(e))

for (assay.name in rownames(new.assays)) {
  assay(e, assay.name) <- reroll(new.assays[assay.name, ], e)
}

# fixup superflouous list for scalar assays
for (u in names(assays(e))) {
  if (class(assay(e, u)[[1]]) %in% c("numeric", "integer")) {
    x <- assay(e, u)
    assay(e, u) <- reroll(unlist(unroll(x)), x)
  }
}

saveRDS(e, file = outfile.path)
#
