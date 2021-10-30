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

e <- readRDS(params$input_file)

cat(sprintf("Expriement file = '%s'\nN genes read = %d\n", params$input_file, nrow(e)))

e <- e[, params$sample]

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

# create null cases by selecting a random location with the gene from a distribution that approximates the observed change points
all_cpts <- unlist(lapply(assay(e,"cpt"), unlist))

random_cp <- Vectorize(function(w) trunc(quantile(all_cpts[all_cpts < w], runif(1))), vectorize.args = "w")
cp_null <- rowRanges(e)[cp$tx_name]
cp_null <- GRanges(seqnames = seqnames(cp_null), 
                   ranges = IRanges(start(cp_null) + random_cp(width(cp_null)), width = 1), 
                   strand = strand(cp_null),
                   tx_name = cp_null$tx_name)
cp_null$is_cp <- 0
cp <- c(cp, cp_null)
cp$is_cp <- as.factor(cp$is_cp)



cp$dna <- Views(BSgenome.Scerevisiae.UCSC.sacCer3::Scerevisiae, resize(cp, width = as.integer(params$outer_window)))
mcols(cp) <- cbind(mcols(cp), alphabetFrequency(cp$dna, as.prob = TRUE)[,(c("A","T"))])
cp <- join_nearest(cp, nuc_centers, distance = TRUE)
cp <- as_tibble(cp)

############################################################
# Run the models
############################################################

h2o.init(nthreads = -1, max_mem_size = "8G", port = 54321)  

# variable names to ignore
v.ignore <- c("seqnames", "start", "end", "width", "strand", "tx_name", "seq_index", "m", "v", "is_cp", "dna", "nuc_id", "gene", "acc", "gene_pos")
x <- setdiff(colnames(cp), v.ignore)
y <- "is_cp"

# Don't pass through the non-atomic columns seqname, strand, dna
# h2o doesn't like it
data <- as.h2o(cp[, - which(colnames(cp) %in% c("seqnames", "strand", "dna"))])
  
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
                              # nfolds = 16,
                              validation_frame = validation,  #only used if stopping_rounds > 0
                              seed = ifelse(params$seed == 0, -1, params$seed))
  
timestamp(sprintf("fit auc = %f,  gini Coef = %f", h2o.auc(rf_fit), h2o.giniCoef(rf_fit)))    
  
  
rf_perf <- h2o.performance(model = rf_fit, newdata = test)
  
# Print model performance
timestamp(sprintf("test auc = %f,  gini Coef = %f", h2o.auc(rf_perf), h2o.giniCoef(rf_perf)))    
  
h2o.shutdown()
sessionInfo()

