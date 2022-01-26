# PickleOneGene.R
# Generage pickle file for one gene

library(rtracklayer)
library(reticulate)
library(GenomicRanges)

gl_filename <- "/n/groups/churchman/rds19/data/S005/genelist.gff"
target_file <- "~/Downloads/target.txt"
gl <- import(gl_filename)
target_gene_name <- "YAR015W"
target_gene <- gl[gl$Name == target_gene_name,]
# hard wired for + strand
bgf <- rtracklayer::import("/n/groups/churchman/rds19/data/S005/fastp-nomito/wt-1.pos.bedgraph.gz")
strand(bgf) <- "+"
ov <- findOverlaps(target_gene, bgf)
scores <- bgf[subjectHits(ov)]
z <- gaps(scores, start = start(target_gene), end = end(target_gene))
z <- z[strand(z) == strand(target_gene) & seqnames(z) == seqnames(target_gene),]
z$score <- 0
result <- sort(c(scores,z))
rscore <- rep.int(result$score, width(result))
writeLines(as.character(rscore), target_file)
