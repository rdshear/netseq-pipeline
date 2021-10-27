# ExploreRF.R
# TODO: DESCRIBE SCRIPT

suppressPackageStartupMessages({
  library(h2o)
  library(optigrab)
  library(SummarizedExperiment)
  library(Biostrings)
  library(GenomicRanges)
  library(rtracklayer)
  library(GEOquery)
  library(tidyverse)
  library(plyranges)
})

rm(list = ls())

default_options <- list(
  input_file = "/n/groups/churchman/rds19/data/S005/HresSE_1024_by_6.rds",
  output_dir = "/n/groups/churchman/rds19/data/S005/h2o/",
  n_genes = 0,
  group = "tmp", # This is the directory under which the results are stored.
  sample = "wt-1",
  feature_mode = "C",
  model_mode = "Kmers/2",
  inner_window = "8",
  outer_window = "16",
  seed = 20190722,
  narrowest_segment = 128
)

if (exists("command_line")) {
  params <- opt_fill(default_options, opts = strsplit(command_line, " ")[[1]])
} else {
  params <- opt_fill(default_options)
}

cat("Parameters\n")
for (i in seq_along(params)) {
  n <- names(params)[[i]]
  cat(n, rep(" ", 20 - nchar(n)), paste(params[[i]], collapse  = ","), "\n", sep = "")
}

if (params$seed != 0) {
  set.seed(params$seed)
}



# given region u and genome g, return array of counts by NT species
get_acgt <- function(u, g) {
  
  Views(g[[as.character(seqnames(u))]], start = start(u), end = end(u)) %>%
    as.character() %>%
    strsplit(., "") %>% 
    table() %>% 
    (function(v) {x <- c(A = 0, C = 0, G = 0, T = 0); x[names(v)] <- v; x}) -> v

  if (as.character(strand(u)) == "-")
    v <- structure(v[4:1], names = names(v))
  v
}

# given gene u and genome g, return GC%
# get_gc <- function(u, g) {
#   v <- get_acgt(u, g)
#   structure(sum(v[c("C", "G")]) / sum(v), names = "PctGC")
# ## TODO USE THIS
todo.fixthis <- function(r) {
  genome_seq <- BSgenome.Scerevisiae.UCSC.sacCer3::Scerevisiae
  r <- split(r, as.character(seqnames(r)))
  v <- mapply(function(v, n) {
    as.character(Views(genome_seq[[n]], start = start(v), end = end(v)))
  }, r, names(r), USE.NAMES = FALSE)
  r <- unlist(r)
  r$ntseq <- unlist(v)
  
  # transform types for Crick (-) strand
  r <- split(r, as.character(strand(r)))
  r$`-`$ntseq <- as.character(reverseComplement(DNAStringSet(r$`-`$ntseq)))
  r <- sort(unlist(r))
  
  names(r) <- paste(r$id, r$seq, r$type, sep = ".")
  r
}




e <- readRDS(params$input_file)

cat(sprintf("Expriement file = '%s'\nN genes read = %d\n", params$input_file, nrow(e)))

# TODO: unroll and allow multiple ... or go ro consensus cp's ?
e <- e[, params$sample]
if (ncol(x) != 1) {
  stop("only one sample allowed")
}

# get rid of genes with no changepoints
e <- e[as.vector(assay(e, "k") > 1),]

# if params$n_genes == 0, then take all the genes, otherwise sample is of size params$n_genes
if (params$n_genes > 0 && params$n_genes < nrow(e)) {
  e <- e[sample(nrow(e), params$n_genes)]
}
cat(sprintf("genes after trimming by --n_genes parameter = %d\n", nrow(e)))



# TODO: MOVE THIS TO POINT OF GENERATION
si <- seqinfo(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene::TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
seqinfo(rowRanges(e)) <- si

Sys.setenv(VROOM_CONNECTION_SIZE = 655360)
supfilename <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE61nnn/GSE61888/suppl/GSE61888_nucs_normed.csv.gz"

q <- read_csv(supfilename, col_names = FALSE, n_max = 2)
names(q) <- as.character(seq_along(q))
q <- q[c(1:7,  which(q[2,] == 0))]
q <- q[, -grep("^Input", q[1,])]
qnames <- unlist(c(q[2, 1:7], q[1,8:ncol(q)]))

q <- read_csv(supfilename, col_names = FALSE, skip = 2, col_select = as.integer(names(qnames)))
colnames(q) <- qnames
nuc_centers <- GRanges(seqnames = seqnames(si)[q$chr], ranges = IRanges(start = q$center, width = 1), seqinfo = si)
mcols(nuc_centers) <- q[, c(-2,-3)]

# Get all the changepoints.
# Drop the first segment and then pull the start position of the rest of them 
cp <- unlist(GRangesList(unlist((assay(e, "segments")))))
cp <- cp[cp$seq_index != 1]
cp <- resize(cp, fix = "start", width = 1)
cp$is_cp <- 1

# create null cases by selecting a random location with the gene
cp_null <- rowRanges(e)[cp$tx_name]
cp_null <- GRanges(seqnames = seqnames(cp_null), 
                   ranges = IRanges(start(cp_null) + sample(width(cp_null), 1), width = 1), 
                   strand = strand(cp_null),
                   tx_name = cp_null$tx_name)
cp_null$is_cp <- 0
cp <- c(cp, cp_null)



cp$dna <- Views(genome_seq, resize(cp, width = as.integer(params$outer_window)))
mcols(cp) <- cbind(mcols(cp), alphabetFrequency(cp$dna, as.prob = TRUE)[,(c("A","T"))])
cp <- join_nearest(cp, nuc_centers, distance = TRUE)

############################################################
# Run the models
############################################################

h2o.init(nthreads = -1, max_mem_size = "8G", port = 54321)  

  # CAUTION: create_roi is memoized 
  feature_ranges <- create_roi(ranges = roi[roi$seq != 1,], 
                               genelist = rowRanges(e), 
                               inner_window = case$inner_window, 
                               outer_window = case$outer_window, 
                               mode = case$feature_mode)

  gr <- cbind(mcols(feature_ranges), kmer2vector(feature_ranges, 1, ""))
  gr$PctGC <- (gr$C + gr$G) / (gr$A + gr$C + gr$G + gr$T)
  gr$PctA <- (gr$A) / (gr$A + gr$C + gr$G + gr$T)
  gr$PctC <- (gr$C) / (gr$A + gr$C + gr$G + gr$T)
  gr$PctG <- (gr$G) / (gr$A + gr$C + gr$G + gr$T)
  gr$PctT <- (gr$T) / (gr$A + gr$C + gr$G + gr$T)
  
  gr$type <- factor(gr$type)
  
  # TODO: Add outer_window
  model_id <- paste(case$model_mode, case$feature_mode, case$inner_window, sep = "_")
  model_id <- sub("/", "-", model_id, fixed = TRUE)
  # response variable name
  y <- "type"
  # variable names to ignore
  v.ignore <- c("seqnames", "start", "end", "width", "strand", "id", "seq", "tss", "stop", "ntseq")

    # Set predictor variables
  x <- character()
  for (mmx in strsplit(case$model_mode, split = "+", fixed = TRUE)[[1]])
  {
    u <- strsplit(mmx,split = "/", fixed = TRUE)[[1]]
    mm <- u[1]
    if (length(u) == 1) {
      mm.param <- 0
    } else {
      mm.param <- as.numeric(u[2])
    }
    x.add <- switch(mm,
      "Kmers" = {
        gr <- cbind(gr, kmer2vector(feature_ranges, mm.param, "Kmer."))
        (function(u) u[grep("^Kmer[.]", u)])(colnames(gr))
        },
      "NcLoc" = {
        # get nucleosome position from Brogaard2012
        # Note most have run liftover_Brogaard.R first
        pathName3 <- "/n/groups/churchman/rds19/data/S003/Brogaard2012/"
        r <- import(paste0(pathName3, "NucleosomeLocations.gff3"))
        t <- findOverlaps(
          resize(shift(feature_ranges, 
                        shift = feature_ranges$cp - start(feature_ranges)), 
                        fix = "center", 
                        width = 1024),
                  r)
        gr$NcLoc <- mapply(function(u, v) {
          if (length(v) == 0) {
            Inf
          } else {
            w <- start(r[v])
            w[which.min(abs(w - u))] - u
          }
        }, feature_ranges$cp, as(t, "list"))
        "NcLoc"
      },
      "NcMod" = {
        psx <- get_histone_mods(feature_ranges)
        gr <- cbind(gr, psx)
        colnames(psx)
      },
      "NcModMax" = {
        psx <- get_histone_mods(feature_ranges, prefix = "NcModMax.", summary_function = max)
        gr <- cbind(gr, psx)
        colnames(psx)
      },
      "PctA" = "PctA",
      "PctC" = "PctC",
      "PctG" = "PctG",
      "PctT" = "PctT",
      "PctGC" = "PctGC",
      "PctACG" =  c("PctA", "PctC", "PctG"),
      stop(sprintf("Invalid model mode value (%s)", mm))
    )
    x <- c(x, x.add)
  }
  data <- as.h2o(gr)
  
  # Partition the data into training, validation and test sets
  splits <- h2o.splitFrame(data = data, 
                           ratios = c(0.7, 0.15),  #partition data into 70%, 15%, 15% chunks
                           seed = 20190725)  #setting a seed will guarantee reproducibility
  train <- splits[[1]]
  validation <- splits[[2]]
  test <- splits[[3]]
  
  timestamp(sprintf("N training rows = %d\nN test rows = %d\n", nrow(train), nrow(test)))
  
  rf_fit <- h2o.randomForest(x = x,
                              y = y,
                              training_frame = train,
                              model_id = model_id,
                              # nfolds = 16,
                              validation_frame = validation,  #only used if stopping_rounds > 0
                              seed = ifelse(params$seed == 0, -1, params$seed))
  
  timestamp(sprintf("fit auc = %f,  gini Coef = %f", h2o.auc(rf_fit), h2o.giniCoef(rf_fit)))    
  
  
  rf_perf <- h2o.performance(model = rf_fit, newdata = test)
  
  # Print model performance
  timestamp(sprintf("test auc = %f,  gini Coef = %f", h2o.auc(rf_perf), h2o.giniCoef(rf_perf)))    
  
  # TODO: create directory for output group if none exists
  saveRDS(list(fit = rf_fit, test = rf_perf, case = case), 
          file = file.path(params$output_dir, params$group,
                           paste0(model_id, ".rds")))
  timestamp(stamp = paste0(model_id, " complete"))
  timestamp()
}

sessionInfo()

