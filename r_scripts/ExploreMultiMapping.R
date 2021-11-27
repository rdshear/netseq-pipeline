# ExploreMultiMapping.R
library(fitdistrplus)
library(SummarizedExperiment)
library(GenomicRanges)
library(GenomicAlignments)
library(magrittr)
library(tidymodels)
library(tidyverse)
library(glue)

e <- readRDS("/n/groups/churchman/rds19/data/S005/HresSE_0_by_5.rds")

variant <- "wt"
e <- e[, which(colData(e)$variant == variant)]
print(
  glue(
    "Variant = {variant}, {ncol(e)} samples selected: {paste0(colnames(e), collapse = ' ')}"
  )
)

zero_run_mask <- function(v, min_width) {
  gr <- range(v)
  w<- gaps(as(v[v$score != 0], "GRanges"), 
         start = start(gr), end = end(gr))
  w <- GenomicRanges::intersect(gr, w)
  w[width(w) >= min_width]
}

na_jaccard <- function(a, b) {
  sum(width(GenomicRanges::intersect(a,b))) / sum(width(GenomicRanges::union(a,b)))
}

na_hamming <- function(a, b) mean(!xor(a, b))

distance_matrix <- function(samples, trial_width) {
  trial_width <- trial_width[1] # TODO: tenp
  x_na <- samples %>% mutate(mask = map(.$score_gpos, zero_run_mask, trial_width))
  distances <- tibble(a = seq(nrow(x_na))) %>% 
    mutate(b = a) %>%
    cross_df(.filter = `>=`) %>%
    mutate(a_mask = map(a, ~x_na$mask[.][[1]]),
           b_mask = map(b, ~x_na$mask[.][[1]]),
          jaccard = map2_dbl(a_mask, b_mask, na_jaccard))
  distances
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
# h <- h[1:8,]
h <- e
g <- as_tibble(rowRanges(h))
range_headers <- colnames(data)
g %<>%
  inner_join(flatten_matrix(h, "scores","score_gpos")) %>%
  inner_join(flatten_matrix(h, "cpt")) %>%
  inner_join(flatten_matrix(h, "mu")) %>%
  inner_join(flatten_matrix(h, "var")) %>%
  inner_join(flatten_matrix(h, "cliff.magnitude")) %>%
  mutate(scores = map(score_gpos, ~as.integer(.$score)),
         scores = map_if(scores, strand == "-", rev),
         cpt = map(cpt, function(u) u[[1]]),
         rle = map(scores, Rle),
         zrl = map(rle, function(z) 
           sort(unique(runLength(z)[runValue(z) == 0]), decreasing = TRUE)),
         max_zrl = map_int(zrl, ~.[1]),
         zero_run_ratio = max_zrl / width) 

x <- g %>% select(Name, sample, cliff.magnitude) %>%
  filter(abs(cliff.magnitude) < 10 & cliff.magnitude !=0) %>%
  pivot_wider(names_from = sample, id_cols = Name, values_from = cliff.magnitude) 

m <- x %>% select(!Name)
  
y <- correlate(m, method = "pearson")
y

xl <- g %>% select(Name, sample, cliff.magnitude) %>%
  filter(abs(cliff.magnitude) < 10 & cliff.magnitude !=0) %>%
  modify_if(is.character, as.factor)

z <- glm(cliff.magnitude ~ ., data = xl)



data_gene <- g %>%
  nest(samples = !c(range_headers)) %>%
  hoist(samples, mask_widths = "zrl", .remove = FALSE) %>% # TODO: .remove for testing only
  # Select candidate mask_widths
  mutate(mask_widths = map(mask_widths, 
                           function(u) {
                             x <- unique(unlist(u))
                             sort(x[x >= minimum_window], decreasing = TRUE)
                           }),
         max_mask_width = map_int(mask_widths, max))
# TODO here's the distance measurement if needed
# data %>%
#   mutate(sim = map2(samples, mask_widths, distance_matrix)) -> data1
# 
#   # look at the detail
#   p <- ggplot(y, aes(x = location, y = score, color = sample)) +
#     geom_jitter(width = 0.25) +
# #    facet_grid(rows = vars(sample)) +
#     labs(title = gene_name) +
#     scale_y_log10()
#   print(p)


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
