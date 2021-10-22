# CombineSummarizedExperiments.R
library(SummarizedExperiment)
library(gplots)
# THis is misnomer. Actually, this is Explore...
# combine two SummarizedExperiments (with concordant rows) into one
# Beginning works. But the ned has specifics around genes and diamide experiments.
combine.SummarizedExperiments <- function(x,y) {
  SummarizedExperiment(assays = mendoapply(function(a, b) cbind(a, b), assays(x), assays(y)),
                       rowRanges = rowRanges(x), 
                       colData = rbind(colData(x), colData(y)))
}

e  <- readRDS("/n/groups/churchman/rds19/data/S005/HresSE_1k_by_6.rds")
# se.screen <- readRDS("/n/groups/churchman/rds19/data/S003/se/exp_screen.rds")
# e <- combine.SummarizedExperiments(se.primary, se.screen)
# e <- se.primary

plot.cpts <- function(object, features = NULL, samples = NULL,
         alpha = 0.4, 
         color = "red", 
         mfrow = NULL,
         label = "", 
         xaxis = TRUE, 
         yaxis = TRUE, ...) {
  oldpar <- par(no.readonly = TRUE)

  e <- object
  if (!is.null(features)) {
    e <- e[features, ]
  }
  if (!is.null(samples)) {
    e <- e[, samples]
  }
  
  par(cex = 0.6, tcl = -0.25, mgp = c(2,0.6, 0),
      mar = c(0.2, 0, 0, 0), oma = c(3, 3, 0.5, 0.5))

  if (is.null(mfrow)) {
    par(mfrow = rev(dim(e)))
  } else {
    par(mfrow = mfrow)
  }
  
  o <- assay(e, "scores")
  p <- assay(e, "changepoints")
  for (s in colnames(e)) {
    for (f in rownames(e)) {
      plot(as.numeric(o[f, s][[1]]), 
           type = "h",
           col = scales::alpha(colour = color, alpha = alpha))
      cp <- c(1, end(p[f, s][[1]]))
      abline(v = cp, col = "blue")
      mtext(paste0(s, " ", f), side = 3, line = -1, adj = 0.015, cex = 0.7, font = 2)
      # TODO: mean by segment, manage axes
    }
  }
  par(oldpar)
  
}

factor.colors <- function(f) {
  f <- as.integer(factor(f))
  n <- max(f)
  n <- n + 3 - n %% 3
  rainbow(n)[as.integer(matrix(seq(n), nrow = 3, byrow = TRUE))[f]]
}

par(mfrow = c(3,1))
col <- factor.colors(colData(e)$variant)

cliff <- assay(e, "cliff.magnitude")
cliff[cliff == 0] <- as.numeric(NA)
m <- assay(e, "mu")

x <- boxplot(cliff, 
        top = TRUE,
        shrink = 0.5,
        las = 2,
        pch = 16, 
        cex = 0.5,
        outline = FALSE,
        col = col,
        notch = TRUE,
        main = "FCP scores by sample", 
        ylab = "First Change Point log fold change")
abline(h = 0, col = "blue")

barplot(colSums(!is.na(cliff)), col = col) 
barplot(colMeans(m), col = col)

odx <- order(apply(cliff, 2, function(u) diff(quantile(u, c(.25, .75), na.rm  = TRUE))))
boxplot(cliff[, odx],
        top = TRUE,
        shrink = 0.5,
        las = 2,
        pch = 16, 
        cex = 0.5,
        outline = FALSE,
        col = col[odx],
        notch = TRUE,
        main = "FCP scores by sample", 
        ylab = "First Change Point log fold change")
abline(h = 0, col = "blue")
barplot(colSums(!is.na(cliff[, odx])), col = col[odx]) 
barplot(colMeans(m[, odx]), col = col[odx])
par(mfrow = c(1,1))


plot(colSums(!is.na(cliff)), colMeans(m))
df1 <- data.frame(m = colMeans(m), N = colSums(!is.na(cliff)))
# Need same sequencing depth to get good answer...maybe we should capture the raw ICL scores?
cor(df1, method = "spearman")

library(ggplot2)
#m1 <- data.frame(assay( e[, c("wt-1", "wt-2", "wt-3")], "fcp"))
m1 <- data.frame(assay( e, "fcp"))
m1 <- m1[which(apply(m1, 1, max) < 800),]
p = ggplot(m1, aes(x = `wt.3`, y = `wt.2`)) +
  ggtitle("FCP for two WT samples") +
  xlab("WT-1 FCP (nt)") +
  ylab("WT-3 FCP (nt)")
p1 = p + 
  geom_point(alpha = 0.5, colour="orange") + 
  geom_density2d() + 
  theme_bw()
# TODO: better contours
plot(p1)


x <- apply(m1, 2, density)
plot(NULL, xlim = c(0, max(sapply(x, function(u) u$x))),
     ylim = c(0, max(sapply(x, function(u) u$y))),
     main = "FCP density")
for (i in seq_along(x)) lines(x[[i]], col = i)
legend("topright", legend = names(x), col = seq_along(x), lty = 1)

# examining YAL005C (SSA1)
# has very long paralog (SSA2, YLL024C)
# but WT-1 .. WT-3 do not show dead spot...why?
# also paralog synonyms may give data on noise level

plot.cpts(e["YAL005C",], mfrow = c(4,6))
x <- assay(e["YAL005C", ], "scores")[1,]
x
lapply(x, function(y) as.integer(y[300:350]))

# No diamide data here
# cd <- colData(e)
# samp <- rownames(cd[cd$group %in% c("primary", "diamide"),])
# plot.cpts(e[c("YAL005C", "YLL024C"), samp], mfrow = c(4,6))

# Let's stick with WT-1 vs DiaT00. IGV / bam confirms that the
# aligned eads just arn't there for DiaT00.
# Check the mean occupancy

x <- assay(e[, c("WT-1", "DiaT00")], "m")
colMeans(x)
colSdDiffs(x)
cor(x, method = "spearman")
plot(x)

wild.samples <-  c("WT-1", "WT-2", "WT-3", "DiaT00")
x <- assay(e[, wild.samples], "m")
pairs(x)
cor(x)

e0 <- e[rowMin(assay(e, "k")[, wild.samples]) > 1, ]

x.cp1 <- assay(e0[, wild.samples], "cp.1")
pairs(x.cp1)
cor(x.cp1)

x <- assay(e0[, wild.samples], "cliff")
x.k <- rowSums(is.finite(x))
table(x.k)

x.k <- rowSums(is.finite(x) & x < 0)
table(x.k)
x.neg <- x[x.k == ncol(x),]
cor(x.neg)
pairs(log2(-x.neg))

x.cp1.reduced <- x.cp1[x.k == ncol(x),]
cor(x.neg)
pairs(log2(-x.neg))

cp1.range <- apply(x.cp1, 1, mad)
hist(cp1.range, breaks = 20)
hist(cp1.range[cp1.range < 100], breaks = 50)
plot(density(cp1.range))
head(table(as.integer(cp1.range)), n = 20)
cp1.matches <- cp1.range[cp1.range < 20]

e1 <- e0[names(cp1.matches), ]
pairs(assay(e1[, wild.samples], "cp.1"))
# BATCH EFFECT!
cor(assay(e1[, wild.samples], "cp.1"))

pairs(assay(e1[, wild.samples], "m"))
# BATCH EFFECT!
cor(assay(e1[, wild.samples], "m"))
boxplot(assay(e1, "cliff"), 
        top = TRUE,
        shrink = 0.5,
        las = 2,
        pch = 16, 
        cex = 0.5,
        outline = FALSE,
        col = col,
        notch = TRUE,
        main = "FCP scores by sample", 
        ylab = "First Change Point log fold change")
abline(h = 0, col = "blue")
plot.cpts(e1[15:25, 1:6])

# TODO: Kill genes with excess zero runs (multimappers)
s <- assay(e1, "scores")
s[] <- sapply(assay(e1, "scores"), function(s) max(runLength(s)))
quantile(as.integer(s), seq(0, 1, .05))
mean(s < 10)
# only 17% of selected genes don't ahve multi-mappers

plot.cpts(e1[15:25, 1:6])
plot.cpts(e["YBR011C", ], mfrow = c(9,8))


# wt mean distributions?
x <- apply(assay(e[, wild.samples], "m"), 2, density)
plot(NULL, xlim = c(0, 5),
     ylim = c(0, max(sapply(x, function(u) u$y))),
     main = "mean occupancy density")
for (i in seq_along(x)) lines(x[[i]], col = i)
legend("topright", legend = names(x), col = seq_along(x), lty = 1)

# Let's look at distance between first 2 cp's

ew <- e[, cd$group %in% c("primary", "screen")]
e.2cp <- sapply(assay(ew, "changepoints"), function(u) {
  ifelse(length(u) < 3, NA, width(u[2]) + 1)
  })
e.2cp <- matrix(e.2cp, nrow = nrow(ew))
rownames(e.2cp) <- rownames(ew)
colnames(e.2cp) <- colnames(ew)
# TODO: make matrix based on means of variants abs(dist) from WT's


y <- apply(e.2cp, 1, function(u) any(is.na(u)))
e.2cp <- e.2cp[!y, ]

heatmap.2(e.2cp, dendrogram = "column",
          na.rm = TRUE,
          ColSideColors = factor.colors(colData(ew)$variant))


rph1.cols <- e.cols$variant %in% c("WT", "rph1")
# Let's try some stats on rph1 vs WT
x <- e.2cp[, rph1.cols]
# eliminate NA's all the way across
y <- apply(x, 1, function(u) all(is.na(u)))
x <- x[!y, ]
