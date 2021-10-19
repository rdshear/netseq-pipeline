# BuildHresSE4Screen.R
# Create HresSE object for NETseq screen data set
# Add derived statistics that are functions of the scores and possibly the changepoints
library(GenomicRanges)
library(rtracklayer)
library(HresSE)
library(fitdistrplus)

set.seed(20190413)

feature.filename <- "/n/groups/churchman/rds19/data/S005/genelist.gff"
occupancy.path <- "/n/groups/churchman/rds19/data/S005/"
changepoint.path <- "/n/groups/churchman/rds19/data/S005/"

outfile.path <- "/n/groups/churchman/rds19/data/S005/hresse/HresSE_all.rds"
sample.group <- "screen"

# create data.frame of samples
# TODO Pull sample table
# u <- c("wt-1", "wt-2")
# y <- lapply(x, function(u) data.frame(variant = u$variant, group = u$group, stringsAsFactors = FALSE))
# z <- Reduce(rbind, y, data.frame())
# z$sample <- names(y)
z <- data.frame(sample = c("wt-1", "wt-2"), variant = c("wt","wt"), group = c("screen", "screen"))
rownames(z) <- z$sample

# TODO: For test only
z <- z[z$group == sample.group, ]

g <- import(feature.filename)
names(g) <- g$ID
# TODO: For test only
g <- sample(g, 32)
 
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
