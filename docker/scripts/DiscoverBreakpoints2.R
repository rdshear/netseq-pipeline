# DiscoverBreakpoints.R
# bedgraph files --> gff3 files
# disregarding multi-mapped regions
# 
# Usage:
# Rscript --vanilla scripts/DiscoverBreakpoints2.R \
#   {config.file.name} \
#   {reference.gene.list.filename} \
#   {max.gene.length} {K.max} {max.genes} \
#   {sample.name} \
#   {output.file}
#
# Example:
# Rscript --vanilla scripts/DiscoverBreakpoints2.R \
    # /n/groups/churchman/rds19/data/S001/refdata/config.json \
    # /n/groups/churchman/rds19/data/S001/refdata/subject_genes.gff3 \
    # 0 \
    # 12 \
    # 5 \
    # wt-1 \
    # /n/groups/churchman/rds19/data/S005/ \
    # /n/groups/churchman/rds19/data/S005/ \

# To run with embedded parameters, set DEBUG.TEST <- TRUE

suppressPackageStartupMessages({
  library(parallel)
#  library(yaml)
  library(GenomicRanges)
  library(rtracklayer)
#  library(SummarizedExperiment)
  library(breakpoint)
  library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
})
set.seed(20190416)

# DEBUG ONLY FROM HERE.....
DEBUG.TEST <- TRUE
if (exists("DEBUG.TEST")) {
  print("DEBUG IS ON -- COMMAND LINE PARAMETERS IGNORED")
  commandArgs <- function(trailingOnly) {
    c("/n/groups/churchman/rds19/data/S001/refdata/config.json",
      "/n/groups/churchman/rds19/data/S001/refdata/subject_genes.gff3",
      "0", # Maximum gene body length
      "12", # Kmax (maximum number of segments)
      "2", # Sample size
      "wt-1",
      "/n/groups/churchman/rds19/data/S005/",
      "/n/groups/churchman/rds19/data/S005/")
  }
}
#  .............TO HERE

args <- commandArgs(trailingOnly = TRUE)

config.filename <- args[1]
subject_genes.filename <- args[2]
GeneMaxLength <- as.numeric(args[3]) # Truncate gene to this length (or inf if 0)
Kmax <- as.numeric(args[4])
maxGenes <- as.numeric(args[5]) # if > 0 sample this number of genes
sample.name <- args[6]
input.directory <- args[7]
output.directory <- args[8]

sprintf("Starting at %s. Sample name = %s. Max genes = %d", 
        Sys.time(), sample.name, maxGenes)

algorithm <- "CEZINB"

options(mc.cores = detectCores())
sprintf("Number of cores detected = %d", getOption("mc.cores"))

infile.pos <- paste0(input.directory, sample.name, ".pos.bedgraph.gz")
infile.neg <- paste0(input.directory, sample.name, ".neg.bedgraph.gz")


sinfo <- seqinfo(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
g <- import(subject_genes.filename)
seqinfo(g) <- sinfo
names(g) <- g$ID

# truncate gene lengths if so desired
if (GeneMaxLength > 0) {
  g <- resize(g, fix = "start", ifelse(width(g) > GeneMaxLength, GeneMaxLength, width(g)))
}

# subset the genes if so desired
if (maxGenes > 0 & maxGenes < length(g)) {
  g <- g[sort(sample(length(g), maxGenes))]
}

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
              source = "DiscoverBreakpoints2",
             algorithm = rep(algorithm, n),
             m = stats["m", ],
             v = stats["v", ])
  result
})

# HACK. Can't set the name in the GRanges constructor

for (i in seq_along(result)) {
  result[[i]]$tx_name = names(result)[i]
}

result <- unlist(GRangesList(result))



#  mcols(result)$tx_name <- as.character(u$tx_name)
  #result$type <- Rle(values = "seq_index", lengths = n)
  
  result
  # # TODO: Carry BIC and logLikelihood
  # # # TODO: get segment statistics
  # # # TODO: reverse negative strand 
  # # IRanges(start = c(1, seg + 1), end = c(seg, length(s)))
  # # # TODO: report multi-map removal areas
  # c(tx_name = u$tx_name, mu = mu, v = v, 
  #   mzl = head(sort(runLength(s), decreasing = TRUE)))

outfile <- file.path(output.directory, paste0(sample.name, "_", algorithm, ".gff3"))
export(result, con = outfile, index = TRUE)

sprintf("Completed at %s\n", Sys.time())
print(sessionInfo())
