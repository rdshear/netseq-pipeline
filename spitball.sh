#!/bin/bash
cd /Users/robertshear/Projects/netseq-pipeline/barcode_tests/umi_tools_from_RX || exit
conda activate umi_tools
time umi_tools dedup -I ../wt-1.aligned.bam --output-stats=deduplicated -S deduplicated.bam \
  --umi-tag=RX --extract-umi-method=tag

time Rscript --vanilla ../../scripts/BamToBedgraph.R . deduplicated.bam wt-1
