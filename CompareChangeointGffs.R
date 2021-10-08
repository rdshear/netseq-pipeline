# Comapre changepoint gffs
library(rtracklayer)
a <- import("~/Downloads/8335801d-447b-4db8-bd23-89efb5a86f5d_netsq_to_changepoint_6a06c5ff-0bb0-4b7a-afd7-22883e5b002e_call-DiscoverBreakpoints_wt-2.cp.CEZINB.gff3")
b <- import("~/Downloads/wt-2.cp.CEZINB.gff3")
g <- import("/n/groups/churchman/rds19/data/S005/genelist.gff")
names(g) <- g$ID

a1 <- split(a, a$tx_name)
b1 <- split(b, b$tx_name)
z <- mapply(identical, a1, b1)
x <- names(z[!z])
# all the genes with differences in changepoints
xg <- g[x]
cat(sprintf("%d regions,  %d (%0.2f%%) differences", length(z), sum(!z), mean(!z)*100))

k_vs_k <- table(sapply(a1, length), sapply(b1, length))
kableExtra::kable(k_vs_k, format = "pipe", )



# xgout <- GRanges(seqnames = seqnames(xg), ranges = IRanges(start = start(xg), end = end(xg)))
# export.bed(xgout, "~/Downloads/fubar.bed")
