# get gene list
library(tidyverse)
# TODO: Description

src = "https://sgd-prod-upload.s3.amazonaws.com/latest/saccharomyces_cerevisiae.gff.gz"
dst = "/n/groups/churchman/rds19/data/S005/genelist.gff"
q <- rtracklayer::readGFF(src, filter = list(type = "gene"))
r <- q[q$orf_classification == "Verified" & q$seqid != "chrmt",]
export.gff3(r, dst)
