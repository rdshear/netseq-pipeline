# fcp_graphic_example.R
# prepare graphic example for proposal comparing tiles, sliding windows and first change points

# TODO: revise to use gff3's or HresSE

rm(list = ls())

library(GenomicFeatures)
library(EBS)
library(parallel)

sample <- "WT-1"


library(GenomicRanges)
library(rtracklayer)

path <- "/n/groups/churchman/rds19/data/"
path1 <- paste0(path, "S001/aligned")
path2 <- paste0(path, "S002/features")

pos.bedgraph <- file.path(path1, paste0(sample, "_pos_cdx.bedgraph.gz"))
neg.bedgraph <- file.path(path1, paste0(sample, "_neg_cdx.bedgraph.gz"))
output.gff3 <- file.path(path2, paste0(sample, "_cliff.gff3"))

bedgraph2rle <- function(infile) {
  p <- mcolAsRleList(import(infile), "score")
  runValue(p) <- ifelse(is.na(runValue(p)), 0, runValue(p))
  # right pad with zeros to chromosome length of reference genome
  sl <- seqlengths(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene::TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
  as(sapply(names(p), function(u) c(p[[u]], Rle(0, sl[u] - sum(runLength(p[u]))))), "RleList")
}

pn <- c(`+` = bedgraph2rle(pos.bedgraph), `-` = bedgraph2rle(neg.bedgraph))
names(pn) <- strand(c("+","-"))

s <- import(output.gff3)
s$mean <- as.numeric(s$mean)
s$variance <- as.numeric(s$variance)
s$fano <- s$variance / s$mean

cds_segment <- s[s$type == "cds_segment"]
cds_ext <- s[s$type == "cds_ext"]


sig.figs <- 5
plotcolors <- c("blue", "red", "darkgreen", "purple")

par(mfcol = c(2, 2), pch = 20)

u <- split(cds_segment, f = cds_segment$seq, drop = FALSE)
lapply(seq_along(u), function(v) {
  x <- u[v][[1]]
  plot(x$mean, x$fano, log = "xy", 
       main = paste0("Sample=", sample, " Segment=", names(u)[v]), 
       xlab = expression(paste("mean occupancy")), 
       ylab = expression(paste("fano factor occupancy")),
       col = "grey")
  abline(a = 0, b = 1, col = "red")
  abline(h = 1, col = "blue")
  plot(x$mean, x$variance, log = "xy", 
       xlab = expression(paste("mean occupancy")), 
       ylab = expression(paste("variance occupancy")),
       col = "grey")
  abline(a = 0, b = 1, col = "red")
  abline(h = 1, col = "blue")
  #  legend("topleft", legend = c("gene", "x=y", "y=1"), lty = c(-1, 1, 1), lwd = 1, pch = c(20, -1, -1), col = c("grey", "red", "blue"))
#  plot((x$variance - x$mean) / x$variance, x$mean^2 / (x$variance - x$mean), 
#       log = "xy", xlab = "p", ylab = "r", col = "grey", pch = 20)
#  abline(a = 0, b = 1, col = "red")
#  plot(x$mean, (x$variance - x$mean) / x$mean^2,
#       log = "xy", xlab = "mean", ylab = expression(paste(alpha, "(dispersion)")), col = "grey", pch = 20)
#  abline(h = 1, col = "blue")
#  abline(a = 0, b = 1, col = "red")
}
)


plot(u$`1`$mean, u$`2`$mean, log = "xy", col = "grey", main = "Segement 1 vs Segment 2", xlab = "mean 1", ylab = "mean 2")
abline(h = 1, col = "blue")
abline(a = 0, b = 1, col = "red")
plot(u$`1`$variance, u$`2`$variance, log = "xy", col = "grey", xlab = "variance 1", ylab = "variance 2")
abline(h = 1, col = "blue")
abline(a = 0, b = 1, col = "red")
plot(u$`1`$fano, u$`2`$fano, log = "xy", col = "grey", xlab = "Fano factor 1", ylab = "Fano factor 2")
abline(h = 1, col = "blue")
abline(a = 0, b = 1, col = "red")



# stop("only do the first part")
# ## DEBUG ONLY: SMALL SAMPLE
# N <- 30
# s <- s[1:(N*3)]
# computed_breakpoints <- width(cds_segment[cds_segment$seq == 1])
# 
# 
# cl <- makeCluster(detectCores(), type = "FORK")
# 
# segs <- parLapply(cl, cds_ext, function(u) {
#   v <- pn[as.character(strand(u))][[1]][seqnames(u)][[1]][start(u):end(u)]
#   EBSegmentation(as.vector(v), model = 3, Kmax = 5)
# })
# 
# 
# break.points <- parLapply(cl, segs, function(seg) {
#   icl <- EBSICL(seg)
#   k <- icl$NbICL
#   result <- sapply(seq(k - 1), function(u) which.max(EBSDistrib(seg, u, k)))
#   if (class(result) != "integer") {
#     result <- integer(0)
#   }
#   c(result, rep(-1, 5 - length(result)))
# })
# stopCluster(cl)
# 
# plot.cp.lines <- function(i, color = "black") {
#   abline(v = break.points[[i]][which(break.points[[i]] > 0)], col = color, lty = 3)
# }
# 
# fp <- function(xdata, s, fcn) {
#   as.numeric(unlist(apply(s, 1, function(u) {
#   rep(fcn(xdata[u[1]:u[2]], na.rm = TRUE), u[2] - u[1] + 1)
#   })))
# }
# 
# plot.it <- function(x, m, v, title = "") {
#   f <- v / m
#   mat <- cbind(m, v, f)
#   matplot(mat, type = "l", lty = 1, col = plotcolors, log = "y")
#   legend("topright", c(title, "mean", "variance", "fano factor"), lty = c(0, 1, 1, 1), col = c("black", plotcolors))
#   points(x, col = "grey", pch = 20)
# }
# 
# for (idx in seq_along(segs)) {
#   gene.name <- cds_ext$Name[idx]
# #  png(paste0("~/Downloads/", gene.name, ".png"), width = 1600, height = 900, units = "px")
#   pdf(paste0("~/Downloads/", gene.name, ".pdf"), width = 32, height = 16)
#   kMax <- max(which(break.points[[idx]] > 0))
# #  par(mfcol = c(2 + 6 + kMax, 1), mai = c(0.1, 1, 0, 0.5), ylog = TRUE)
#   
#   layout(matrix(c(1, 2, 3, 4, 5, 7, 9, 11, 6, 8, 10, 12), nrow = 4, ncol = 3, byrow = FALSE))
#   par(mai = c(0, 0.3, 0, 0), oma =c(5, 5, 5, 0), ylog = TRUE)
#   
#   x <- segs[idx][[1]]@data
#   plot(x, log = "y", pch = 20, ylab = "occupancy", col = "grey")
#   legend("topright", c("occupancy", "change point"), lty = c(0,3), pch = c(20, -1), col = c("grey", "black"))
#   plot.cp.lines(idx)
#   mtext("NT downstream of TSS", outer = TRUE, line = 2, side = 1)
#   mtext("count", outer = TRUE, line = 2, side = 2)
#   mtext(gene.name, outer = TRUE, line = 1, side = 3)
#   
#   s <- (function(v) {
#     z <- length(v)
#     cbind(v[1:(z - 1)], v[2:z] - 1)})(c(1, (function(u) u[u > 0])(break.points[idx][[1]]), length(x)))
#   #TODO: linear regression on segments here
#   
#   m <- fp(x, s, mean)
#   v <- fp(x, s, var)
#   plot.it(x, m , v, "exact baysiean segmentation")
#   
#   # TODO: TEMP HACK
#   kMax <- 2
#   for (k in 1:kMax) {
#     plot(EBSDistrib(segs[idx][[1]], k,5), type = "l", col = plotcolors[k])
#     legend("topright", paste0("posterior probability segment ", k), col = plotcolors[k], lty = 1)
#     plot.cp.lines(idx)
#   }
# 
#   for (window in 2^seq(5, 8)) {
#     m <- filter(x, rep(1 / window, window))
#     v <- filter(x^2, rep(1 / (window - 1), window))
#     plot.it(x, m, v, paste0("sliding window ", window, " nt"))
#     
#     s <- t(sapply(seq(from = 1, to = length(x) - 1, by = window), function(i) c(i, min(length(x), i + window - 1))))
#     m <- fp(x, s, mean)
#     v <- fp(x, s, var)
#     plot.it(x, m, v, paste0("tiled, width = ", window, " nt"))
#   }
#   dev.off()
# }
# 
