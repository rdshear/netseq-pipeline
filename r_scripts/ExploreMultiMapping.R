# ExploreMultiMapping.R
library(SummarizedExperiment)
library(GenomicRanges)
library(GenomicAlignments)

e <- readRDS("/n/groups/churchman/rds19/data/S005/HresSE_0_by_5.rds")

g <- c("YBR204C")
s <- c("wt-1")

h <- e[g,s]
sg <- assay(h, "cpt")[[1]][[1]]
score <- assay(h, "scores")[[1]]
data <- score$score
if (as.character(strand(rowRanges(h))) == "-") {
  data <- rev(data)  
}

cpt_np <- cpt.np(data, minseglen = 16)
plot(cpt_np)
abline(v = sg, col = "purple")

cpt_meanvar <- cpt.meanvar(data,  method = "PELT", minseglen = 16, penalty = "CROPS", pen.value = c(5,5000))
cpt_meanvar
pen.value.full(cpt_meanvar)
plot(cpt_meanvar, ncpts = 9)


cpt_mean <- cpt.mean(data, method = "PELT", minseglen = 16)
cpt_mean


# looking at distributions and censoring
z <- Rle(data)
y <- runLength(z)[runValue(z) == 0]
hist(y)

mean(data, na.rm = TRUE)
censor.lim <- 6
c.slots <- which(runLength(z) > censor.lim & runValue(z) == 0)
runValue(z)[c.slots] <- NA
cm <-mean(z, na.rm = TRUE)
cm
dpois(0:12, cm)
plot(dpois(0:12, cm))


bamFileName <- "/Users/robertshear/temp_save/cromwell-executions/NETseq/fcdadec9-dd6a-4111-839a-bca83daba896/call-StarAlign/execution/wt-4.aligned.bam"


filter <- rowRanges(h)
reads <- readGAlignments(bamFileName, use.names = TRUE,
                              param = ScanBamParam(which = filter, tag = c("NH","HI","AS", "RX"), what=c("qual", "mapq")))


# bamFileName <- "/Users/robertshear/temp_save/cromwell-executions/NETseq/fcdadec9-dd6a-4111-839a-bca83daba896/call-StarAlign/execution/short.bam"
# reads <- readGAlignments(bamFileName, use.names = TRUE,
#                               param = ScanBamParam(what=c("qual", "mapq"), tag = c("NH","HI","AS", "RX")))
# reads <- reads[names(reads) == 'NB500917:93:HHG5KBGXX:1:23302:13142:17955']

table(mcols(reads)[,c("NH","HI")])


idx <- which(mcols(reads)$HI == 9)
r1 <- reads[idx]
r1


r2 <- readGAlignments("/n/groups/churchman/rds19/data/S005/HresSE_1024_by_6.rds", use.names = TRUE,
                         param = ScanBamParam(what = c(qname = names(r1)), tag = c("NH","HI","AS", "RX")))


r3 <- readGAlignments("/n/groups/churchman/rds19/data/S005/HresSE_1024_by_6.rds", use.names = TRUE,
                      param = ScanBamParam(tag = c("NH","HI","AS", "RX")))
