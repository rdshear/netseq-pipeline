# ViewOrfsImage.R

# Call:
# RScript --vanilla ViewOrfsImage.R \
#   display.window.size 
#   minimum.fcp.location 
#   sample.name 
#   fcp-rdata-file.RData
#   orf-rdata-file.RData
#   output-file.pdf

set.seed(6272017)

library(fields, quietly = TRUE)
library(GenomicRanges, quietly = TRUE)

lower.trim.parameter <- 0.15 # lower 'trim' cutoff, mean occupancy / nt
running.mean.window <- 5  # running mean width for smothing img heat map

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop("Must have 6 command line parameters")
}
display.window <- as.numeric(args[1])
if (is.na(display.window)) {
  stop("display.window.size parameter must be numeric")
}
minimum.fcp.location <- as.numeric(args[2])
if (is.na(minimum.fcp.location)) {
  stop("minimum.fcp.location parameter must be numeric")
}
sample.name <- args[3]
fcp.RData <- args[4]
orf.RData <- args[5]
output.filename <- args[6]


load(fcp.RData)
load(orf.RData)

# TODO: move transpose and -1 --> NA to generating program, also move lower 10% trim to generating program
# ALSO transform -1 --> NA
break.points <- t(break.points)
# remove low expression levels ~ 10% of WT
x <- which(sapply(orf, function(u) mean(u$score[[1]])) < lower.trim.parameter)
orf <- orf[-x,]
break.points <- break.points[-x,]
minus.to.na <- function(x) ifelse(x < 0, NA, x)
break.points <- matrix(minus.to.na(as.vector(break.points)), dim(break.points))
row.names(break.points) <- names(orf)

idx <- apply(break.points[,2:6], 1, function(u) min(u[u > minimum.fcp.location], na.rm = TRUE))

# sort by first break point (FBP) position
# eliminate FBP's that exceed the window size - minimum.fcp.location
idx <- idx[-which((idx + minimum.fcp.location > 
                     sapply(seq_along(orf), function(i) length(orf[i]$score[[1]])))
                 | (idx < minimum.fcp.location))]

idx <- sort(idx[idx < display.window])

img <- sapply(names(idx), function(gene.name) {
  u <- as.numeric(runmean(orf[gene.name]$score[[1]], running.mean.window))
  length(u) <- display.window
  u
})

ymax <- length(img)
xmax <- display.window
# normalize the heat map
img <- scale(img)
# display relative density as heat map
image.plot(x = 1:xmax , y = 1:length(idx), z = log2(img), zlim = c(-4,4), col = rainbow(20), 
      xlab = "NT from TSS", ylab = "gene", legend.lab = "normalized log2 occupancy",
      main = paste("RNAPII Occupancy. Sample", sample.name), 
      sub = "Sorted by position of first change point")
lines(idx, 1:length(idx), col = "black")

# stretch / shrink the tracks so that the CP is in the middle
# mapping image to stretch/shrink

# not very good or very elegant dilation, but okay for preliminary work
img.xformed <- sapply(seq_along(idx), function(gene.id) {
  cp <- idx[gene.id]
  # stretch the first segment
  z <- round(((0:(cp - 1) * display.window / 2 / cp)),digits = 0)
  seg1 <- rep(img[1:cp, gene.id], c(z[2:length(z)], round(display.window / 2,0)) - z)
  # grab the second segment
  seg2.data <- img[(cp + 1):display.window, gene.id]
  seg2.data <- seg2.data[!is.na(seg2.data)]
  # now shrink it
  seg2.idx <- findInterval(seq(length(seg2.data)) * display.window / 2 / length(seg2.data), seq(display.window / 2))
  seg2 <- sapply(1:(display.window / 2), function(u) mean(seg2.data[seg2.idx == u]))
  c(seg1, seg2)
})
image.plot(x = ((1:nrow(img.xformed)) - nrow(img.xformed) / 2) / nrow(img.xformed) * 2, y = 1:ncol(img.xformed), z = log2(img.xformed), 
      zlim = c(-4,4), col = rainbow(20),
      xlab = "Change point shifted to middle", ylab = "gene", legend.lab = "normalized log2 occupancy",
      main = paste("RNAPII Occupancy. Sample", sample.name), 
      sub = "First segment stretched, second segment compressed")
abline(v = 0, col = "black")

# what do the slopes of the first segment look like?
slope.first <- sapply(100:1000, function(i) {
  lo <- 1
  hi <- idx[i] - 1
  coef(glm(x ~ t, data = data.frame(t = lo:hi, x = img[lo:hi,i])))["t"]
  }
)

slope.next100 <- sapply(100:1000, function(i) {
  lo <- idx[i] + 1
  hi <- idx[i] + 100
  coef(glm(x ~ t, data = data.frame(t = lo:hi, x = img[lo:hi, i])))["t"]
  }
)

boxplot(slope.first, slope.next100, names = c("1:FCP-1","FCP+1:FCP+100"),
        main = paste("Intrasegment slopes. Sample", sample.name), 
        sub = "Slope of linear reg within segment", ylab = "slope",
        outline = FALSE, notch = TRUE)
abline(h = 0, col = "blue")

meta.cp.min <- 100
x <- sapply(which(idx > meta.cp.min & idx + meta.cp.min < display.window), function(u) {
  cp <- idx[u]
  img[(cp - meta.cp.min):(cp + meta.cp.min), u]
})
plot(seq(from = -meta.cp.min, to = meta.cp.min), rowMeans(x, na.rm = TRUE), type = 'l',
     main = paste("Mean occupancy aligned to first change point", sample.name),
     xlab = "NT from FCP",
     ylab = "Mean occupancy")

# compute relative expression: before & after
expr.levels <- sapply(names(idx), function(i) {
  u <- orf[i]$score[[1]]
  v <- idx[i]
  c(left = mean(u[1:(v - 2)]), right = mean(u[(v + 2):min(display.window, length(u))], na.rm = TRUE))
  })
rel.expr <- apply(expr.levels, 2, function(u) u[1] / u[2])

plot(log2(rel.expr), idx)
abline(v = 0, col = "red")

plot(log2(rel.expr), log2(expr.levels[1,]))
abline(v = 0, col = "red")

z <- as.integer(orf["YLL024C"]$score[[1]])
plot(z)

hist(log2(rel.expr))
quantile(rel.expr, seq(0,1,.05))
sum(rel.expr > 1) / length(rel.expr)
