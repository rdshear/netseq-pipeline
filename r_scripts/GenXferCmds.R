library(tidyverse)
indat <- read.csv("~/Downloads/fsort.csv", header = FALSE)
colnames(indat) <- c("bytes", "date", "inref")

indat$filename <- sapply(indat$inref, function(u) {
  x <- strsplit(u, "/")[[1]]
  x[length(x)]
}, USE.NAMES = FALSE)

indat <- indat[indat$bytes > 0,]

suf <- c(".Log.final.out",".Log.out",".Log.std.out",".bam",".dedup.log",".neg.bedgraph.gz",".pos.bedgraph.gz")
pat <- paste(suf, collapse = "|")

target <- "gs://fc-secure-eea6de1f-1832-463e-bf86-820192da1f92/20211231-commit243d3d562/"

result <- indat[str_detect(indat$filename, pat),]


result1 <- result[1,]

cmds <- paste0("gsutil cp ", result$inref, " ", target, result$filename)

writeLines(cmds, "~/Downloads/xfers.sh")

indat$jobref <- sapply(strsplit(indat$inref, split = "/"), function(u) u[4])
x <- split(indat$date, indat$jobref)
y <- sort(sapply(x, max))
z <- data.frame(jobroot = paste0("gs://fc-secure-eea6de1f-1832-463e-bf86-820192da1f92/", names(y)), jobdate = y, stringsAsFactors = FALSE)

killdirs <- paste0("gsutil rm -r \"", z$jobroot, "/*\"")
writeLines(killdirs, "~/Downloads/killdirs.sh")
