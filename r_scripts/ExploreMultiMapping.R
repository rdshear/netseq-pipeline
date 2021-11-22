# ExploreMultiMapping.R
library(changepoint.np)
library(fitdistrplus)
library(SummarizedExperiment)
library(GenomicRanges)
library(GenomicAlignments)
library(tidygenomics)
library(magrittr)
library(tidyverse)
library(glue)

e <- readRDS("/n/groups/churchman/rds19/data/S005/HresSE_1024_by_6.rds")

variant <- "wt"
e <- e[, which(colData(e)$variant == variant)]
print(
  glue(
    "Variant = {variant}, {ncol(e)} samples selected: {paste0(colnames(e), collapse = ' ')}"
  )
)

zero_run_mask <- function(s, min_width) {
  x <- Rle(s)
  runValue(x)[which(runLength(x) >= min_width & runValue(x) == 0L)] <- as.integer(NA)
  x <- !is.na(x)
  as(x, "vector")
}

na_hamming <- function(a, b) mean(!xor(a, b))

distance_matrix <-  function(samples, trial_width) {
  x_na <- samples %>% mutate(mask = map(.$scores, zero_run_mask, trial_width))
  m_size <- nrow(x_na)
  m <- matrix(NA, ncol = m_size, nrow = m_size)
  m_idx <- which(upper.tri(m, diag = FALSE), arr.ind=TRUE) 
  for (md in 1:nrow(m_idx)) {
    row <- m_idx[md, 1]
    col <- m_idx[md, 2]
    browser()
    m[row ,col] <- na_hamming(x_na$mask[row][[1]], x_na$mask[col][[1]])
  }
 answer <- tibble_row(gene = gene_name, width = trial_width,
                      score = log10(1 - mean(m, na.rm = TRUE)))
  print(trial_width)
  print(m)
  print(log10(1 - mean(m, na.rm = TRUE)))
  answer
}

minimum_window <- 5

mu <- assay(e, "mu")

density.from <- 1
density.to <- 5
# select moderately dense genes
# for simiplicty the samples will all have to be within the range

# could use tidyExperiment if they fix the invisible cast
flatten_matrix <- function(experiment, assay_name, column_name = assay_name) {
  assay(experiment, assay_name) %>%
    as_tibble(rownames = "Name") %>%
    pivot_longer(cols = !Name, names_to = "sample", values_to = column_name)
}

h <- e[as.vector(rowMin(mu) > density.from & rowMax(mu) < density.to),]
# TODO DEBUG ONLY
h <- h[1:8,]

data <- as_tibble(rowRanges(h))
range_headers <- colnames(data)
data %<>%
  inner_join(flatten_matrix(h, "scores","score_gpos")) %>%
  inner_join(flatten_matrix(h, "cpt")) %>%
  mutate(scores = map(score_gpos, ~as.integer(.$score)),
         scores = map_if(scores, strand == "-", rev),
         cpt = map(cpt, function(u) u[[1]]),
         rle = map(scores, Rle),
         zrl = map(rle, function(z) 
           sort(unique(runLength(z)[runValue(z) == 0]))),
         max_zrl = map_int(zrl, ~.[1]),
         zero_run_ratio = max_zrl / width) %>%
  nest(samples = !c(range_headers)) %>%
  hoist(samples, mask_widths = "zrl", .remove = FALSE) %>% # TODO: .remove for testing only
  # Select candidate mask_widths
  mutate(mask_widths = map(mask_widths, 
                           function(u) {
                             x <- unique(unlist(u))
                             sort(x[x >= minimum_window])
                           }))



final_answer <- lapply(seq(nrow(scores)), function(i) {
  gene_name <- rownames(scores)[i]
  candidate_widths <- sort(unique(unlist(map(scores[i,], function(u) u$zrl))))
  candidate_widths <- candidate_widths[candidate_widths > minimum_window]
})




hist(log10(x$z_density))
# look at very long zero-lits
top_zrl <- sapply(scores, function(u)
  u$zrl[1])
dim(top_zrl) <- dim(h)
dimnames(top_zrl) <- dimnames(h)

very_wide_zero_runs <- which(rowMax(top_zrl) > 100)
top_zrl[very_wide_zero_runs, ]

for (i in very_wide_zero_runs) {
  gene_name <- rownames(scores)[i]
  
  x <- sapply(scores[i, ], function(u)
    u$score)
  colnames(x) <- colnames(scores)
  y <- as_tibble(x) %>%
    mutate(location = seq(nrow(.))) %>%
    pivot_longer(cols = starts_with("wt-"),
                 names_to = "sample",
                 values_to = "score")
  
  # look at the detail
  p <- ggplot(y, aes(x = location, y = score, color = sample)) +
    geom_jitter(width = 0.25) +
#    facet_grid(rows = vars(sample)) +
    labs(title = gene_name) +
    scale_y_log10()
  print(p)
}

# TODO: FORCE TO DATA FRAME
# cpt_np <- cpt.np(data, minseglen = 16)
# plot(cpt_np)
# abline(v = sg, col = "purple", lwd = 3)
#
# cpt_meanvar <- cpt.meanvar(data,  method = "PELT", minseglen = 16, penalty = "CROPS", pen.value = c(5,5000))
# cpt_meanvar
# pen.value.full(cpt_meanvar)
# plot(cpt_meanvar, ncpts = 9)
#
#
# cpt_mean <- cpt.mean(data, method = "PELT", minseglen = 16)
# cpt_mean
#
#
# # looking at distributions and censoring
# z <- Rle(data)
# y <- runLength(z)[runValue(z) == 0]
#
# mean(data, na.rm = TRUE)
# censor.lim <- 6
# c.slots <- which(runLength(z) > censor.lim & runValue(z) == 0)
# runValue(z)[c.slots] <- NA
# cm <-mean(z, na.rm = TRUE)
# cm
# dpois(0:12, cm)
# plot(dpois(0:12, cm))
#
#
# bamFileName <- "/Users/robertshear/temp_wt-4/outputs/wt-4.dedup.bam"
#
#
# filter <- rowRanges(h)
# reads <- readGAlignments(bamFileName, use.names = TRUE,
#                               param = ScanBamParam(which = filter, tag = c("NH","HI","AS", "RX"), what=c("qual", "mapq")))
#

# bamFileName <- "/Users/robertshear/temp_save/cromwell-executions/NETseq/fcdadec9-dd6a-4111-839a-bca83daba896/call-StarAlign/execution/short.bam"
# reads <- readGAlignments(bamFileName, use.names = TRUE,
#                               param = ScanBamParam(what=c("qual", "mapq"), tag = c("NH","HI","AS", "RX")))
# reads <- reads[names(reads) == 'NB500917:93:HHG5KBGXX:1:23302:13142:17955']
#
# table(mcols(reads)[,c("NH","HI")])
#
#
# idx <- which(mcols(reads)$HI == 9)
# r1 <- reads[idx]
# r1
#
#
# r2 <- readGAlignments("/n/groups/churchman/rds19/data/S005/HresSE_1024_by_6.rds", use.names = TRUE,
#                          param = ScanBamParam(what = c(qname = names(r1)), tag = c("NH","HI","AS", "RX")))
#
#
# r3 <- readGAlignments("/n/groups/churchman/rds19/data/S005/HresSE_1024_by_6.rds", use.names = TRUE,
#                       param = ScanBamParam(tag = c("NH","HI","AS", "RX")))
