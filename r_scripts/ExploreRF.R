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
  model_name = "cpa_permute",
  input_file = "/n/groups/churchman/rds19/data/S005/HresSE_0_by_6.rds",
  output_dir = "/n/groups/churchman/rds19/data/S005/h2o/",
  n_genes = 0,
  group = "tmp", # This is the directory under which the results are stored.
  sample = c("wt-1", "wt-2", "wt-3", "wt-4"),
  inner_window = "8",
  outer_window = "16",
  seed = 20190722,
  null_case_method = "permute" # freq for cp location ~ observed cp offset from TSS, permute for permuting modification levels within nucleosome postion
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

# the gene information in this file appears to be unreliable. We will construct our own nucleosome posotion information
#### WIP #####
r <- rowRanges(e)
x <- findOverlaps(r, nuc_centers) %>%
  as_tibble() %>%
  mutate(center = start(nuc_centers[subjectHits]), 
         strand = as.character(strand(r[queryHits]))) %>%
  mutate(sortorder = if_else(strand == "-", -center, center)) %>%
  group_by(queryHits) %>%
  arrange(sortorder, .by_group = TRUE) %>%
  mutate(gene_pos = row_number()) %>%
  ungroup()


nuc_centers <- nuc_centers[x$subjectHits]
nuc_centers$gene_pos <- x$gene_pos


# Get all the changepoints ... unroll the SE matrix
# omit genes without changepoints and add sample name and variant to the list
cp <- do.call(c, lapply(seq(ncol(e)), function (u) {
        elcl <- e[,u]
        x <- unlist(GRangesList(unlist(assay(elcl, "segments"))[as.vector(assay(elcl, "k") > 1)]))
        x$sample <- colnames(elcl)
        x$variant <- colData(elcl)$variant
        x
}))

# Drop the first segment and then pull the start position of the rest of them 
cp <- cp[cp$seq_index != 1]
cp <- resize(cp, fix = "start", width = 1)
cp$is_cp <- 1

switch (params$null_case_method,
  # create null cases by selecting a random location with the gene from a distribution that approximates the observed change points
  "freq" = {
    all_cpts <- unlist(lapply(assay(e,"cpt"), unlist))
  
    random_cp <- Vectorize(function(w) trunc(quantile(all_cpts[all_cpts < w], runif(1))), vectorize.args = "w")
    cp_null <- rowRanges(e)[cp$tx_name]
    cp_null <- GRanges(seqnames = seqnames(cp_null),
                       ranges = IRanges(start(cp_null) + random_cp(width(cp_null)), width = 1),
                       strand = strand(cp_null))
    mcols(cp_null) <- mcols(cp)
    cp_null$is_cp <- 0
    nuc_centers_null <- nuc_centers
  },

# create the null case by permuting the nuc_centers values, but only within nucleosome position (to avoid distsnce from TSS bias)
  "permute" = {
    cp_null <- cp
    cp_null$is_cp = 0
    # TODO Tier by position
    nuc_centers_null <- nuc_centers
    mcols(nuc_centers_null) <- mcols(nuc_centers)[sample(length(nuc_centers), length(nuc_centers)),]

    },
#
  stop("null case menthod not found")
)

cp <- c(join_nearest_downstream(cp, nuc_centers, distance = TRUE),
        join_nearest_downstream(cp_null, nuc_centers_null, distance = TRUE))
cp$dna <- Views(BSgenome.Scerevisiae.UCSC.sacCer3::Scerevisiae, resize(cp, width = as.integer(params$outer_window)))
mcols(cp) <- cbind(mcols(cp), alphabetFrequency(cp$dna, as.prob = TRUE)[,(c("A","T"))])
cp$is_cp <- as.factor(cp$is_cp)

cp$dna <- as.character(cp$dna)
cp <- as_tibble(cp)

############################################################
# Run the models
############################################################

h2o.init(nthreads = -1, max_mem_size = "8G", port = 54321)  

# variable names to ignore
v.ignore <- c("seqnames", "start", "end", "width", "strand", "tx_name", 
              "seq_index", "m", "v", "is_cp", "dna", "nuc_id", "gene", 
              "acc", "gene_pos")
x <- setdiff(colnames(cp), v.ignore)
y <- "is_cp"

data <- as.h2o(cp, destination_frame = params$model_name)
  
# Partition the data into training, validation and test sets
splits <- h2o.splitFrame(data = data, 
             ratios = c(0.7, 0.15),  #partition data into 70%, 15%, 15% chunks
             destination_frames = 
               paste0(params$model_name, "_",  c("fit", "validate", "test")),
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
  
if (!interactive()) {
  h2o.shutdown(prompt = TRUE)
}

sessionInfo()

