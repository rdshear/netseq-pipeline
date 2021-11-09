# DiscoverMotifs.R

suppressPackageStartupMessages({
  library(optigrab)
  library(SummarizedExperiment)
  library(BSgenome.Scerevisiae.UCSC.sacCer3)
  library(rGADEM)
  library(TFBSTools)
  library(JASPAR2018)
  library(tidyverse)
  library(plyranges)
})


rm(list = ls())

default_options <- list(
  input_file = "/n/groups/churchman/rds19/data/S005/HresSE_0_by_6.rds",
  output_file = "/n/groups/churchman/rds19/data/S005/motifs",
  n_genes = 512,
  sample = c("wt-1", "wt-2", "wt-3", "wt-4"),
  seed = 20211006
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

# TODO Pick the highest quality cp's
# will this do well enough?
k <- assay(e, "k")
x <- which(rowMax(k) == rowMin(k) & (rowMax(k) > 1) & rowMin(assay(e, "mu")) > 0.3)
e <- e[x,]

# if params$n_genes == 0, then take all the genes, otherwise sample is of size params$n_genes
if (params$n_genes > 0 && params$n_genes < nrow(e)) {
  e <- e[sample(nrow(e), params$n_genes)]
}
cat(sprintf("genes after trimming by --n_genes parameter = %d\n", nrow(e)))

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
cp <- flank(cp, width = 20, both = TRUE)
strand(cp) <- "*"
seqs <- Views(BSgenome.Scerevisiae.UCSC.sacCer3::Scerevisiae, cp)
seqs <- as(seqs, "XStringSet")
# get the sequences


#seqs <- Views(Scerevisiae, cp)

gadem <- GADEM(seqs, genome = Scerevisiae, verbose = TRUE,
           numTop3mer = 12, numTop4mer = 36, numTop5mer = 42)
novel_motifs <- consensus(gadem)
plot(gadem[1])
consensusMotifs <- consensus(gadem)
print(consensusMotifs)

unknown_motif <- getPWM(gadem)[[1]]

unknown_pwm <- PWMatrix(profileMatrix = unknown_motif)

pwm_library <- getMatrixSet(JASPAR2018, opts = list(collection = "CORE",
                      species = "Saccharomyces cerevisiae", matrixtype = "PWM"))

pwm_sim = PWMSimilarity(pwm_library, unknown_pwm, method = 'Pearson')

pwm_library_list = lapply(pwm_library, function(x){
  data.frame(ID = ID(x), name = name(x))
})

# combine the list into one data frame
pwm_library_dt = dplyr::bind_rows(pwm_library_list)

# fetch the similarity of each motif to our unknown motif
pwm_library_dt$similarity = pwm_sim[pwm_library_dt$ID]

# find the most similar motif in the library
pwm_library_dt = pwm_library_dt[order(-pwm_library_dt$similarity),]

head(pwm_library_dt)
