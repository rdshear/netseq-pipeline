# ExploreAutoML4.R

suppressPackageStartupMessages({
  library(HresSE)
  library(Biostrings)
  library(GenomicRanges)
  library(rtracklayer)
  library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
  library(BSgenome.Scerevisiae.UCSC.sacCer3)
  library(tidyverse)
  library(h2o)
  library(optigrab)
})

rm(list = ls())

gr <- import( "/n/groups/churchman/rds19/data/S004/mark_features/screen//peakfeature_WT-1.gff3") 
gr %>% transform(dwell = as.integer(dwell), seq = as.integer(seq), 
            score = as.numeric(score), start = as.integer(start), 
            end = as.integer(end), width = as.integer(width), 
            is_cp = as.logical(is_cp)) %>%
  as_tibble %>% print(width = Inf)  -> w
# %>%
#   dplyr::select(ID, seq, seqnames, start, strand, is_cp, mark, dwell, score) %>%
#   print(width = Inf) %>%
#   pivot_wider(names_from = mark, values_from = c(dwell, score)) -> w



w_real <- filter(w, is_cp)

w_real %>% group_by(mark) %>% 
  ggplot(aes(x = dwell, group = mark, title = mark)) + facet_wrap(. ~ mark) +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.9) +
    geom_vline(xintercept = c(-55,0, 55))



g <- BSgenome.Scerevisiae.UCSC.sacCer3
gr <- import( "/n/groups/churchman/rds19/data/S004/mark_features/screen//peakfeature_WT-1.gff3")
grx <- split(gr, gr$is_cp)
for (w in 1:10) {
  at <- letterFrequency(getSeq(g, resize(grx$TRUE, w)), letters = c("A", "T"), as.prob = TRUE)
  cat(w, "  ", colMeans(at), "\n")
}

gs <- function(gr1, w, n, prefix) {
  u <- getSeq(g, resize(gr1, w)[sample(length(gr1), n)])
  u <- strsplit(as.character(u), "")
  names(u) <- paste0(prezfix, seq_along(u))
  u <- t(as_tibble(u))
  colnames(u) <- paste0("V", seq(ncol(u)))
  as_tibble(u)
}

window <- 13
N <- 500
background <- gs(grx$`FALSE`, window, N, "BG")
target <- gs(grx$`TRUE`, window, N, "FG")
indata <- rbind(cbind(background, response = rep(0, nrow(background))),
  cbind(target, response = rep(1, nrow(target))))

library(h2o)
h2o.init()
