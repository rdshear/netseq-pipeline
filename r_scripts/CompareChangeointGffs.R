# Comapre changepoint gffs

# TODO Revisit at end to make table of thesis vs paper differences
library(rtracklayer)
a <- import("/n/groups/churchman/rds19/data/S004/se/WT-2_CEZINB.gff3.bgz")
b <- import("~/Downloads/wt-2.cp.CEZINB.gff3")
g <- import("/n/groups/churchman/rds19/data/S005/genelist.gff")
names(g) <- g$ID

a1 <- split(a, a$tx_name)
b1 <- split(b, b$tx_name)

shared.names <- intersect(names(a1), names(b1))
a1 <- a1[shared.names]
b1 <- b1[shared.names]

z <- mapply(identical, a1, b1)
x <- names(z[!z])
# all the genes with differences in changepoints
xg <- g[x]
cat(sprintf("%d regions,  %d (%0.2f%%) differences", length(z), sum(!z), mean(!z)*100))

k_vs_k <- table(sapply(a1, length), sapply(b1, length))
kableExtra::kable(k_vs_k, format = "pipe", )



# xgout <- GRanges(seqnames = seqnames(xg), ranges = IRanges(start = start(xg), end = end(xg)))
# export.bed(xgout, "~/Downloads/fubar.bed")
