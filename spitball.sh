#!/bin/bash
cd /Users/robertshear/Projects/netseq-pipeline/barcode_tests/umi_tools_from_RX || exit
conda activate umi_tools
time umi_tools dedup -I ../wt-1.aligned.bam --output-stats=deduplicated -S deduplicated.bam \
  --umi-tag=RX --extract-umi-method=tag

time Rscript --vanilla ../../scripts/BamToBedgraph.R . deduplicated.bam wt-1

gzcat /n/groups/churchman/GSE159603/wt-1.fastq.gz \
| perl "$(which prinseq-lite.pl)" -fastq stdin \
        -out_good stdout \
		-out_bad output.bad \
		-no_qual_header \
		-min_len 7 -min_qual_mean 20 -trim_right 1 -trim_ns_right 1 \
		-trim_qual_right 20 -trim_qual_type min \
		-trim_qual_window 1 -trim_qual_step 1 \
		2>> log


gzcat cromwell-executions/NETseq/182e220f-837c-4b35-930d-680c7c819426/call-StarAlign/inputs/-1850460386/xwt-1.fastq.gz \
| head -n 400 | perl "$(which prinseq-lite.pl)" -fastq stdin \
		-out_bad output.bad \
		-no_qual_header \
		-min_len 7 -min_qual_mean 20 -trim_right 1 -trim_ns_right 1 \
		-trim_qual_right 20 -trim_qual_type min \
		-trim_qual_window 1 -trim_qual_step 1 \
		2>> log