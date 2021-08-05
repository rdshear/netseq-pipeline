#! /bin/bash
cd ~/temp || exit
conda activate cpa

cromwell run -i ~/Projects/netseq-pipeline/test/bigtest_inputs.json -t wdl \
    -o /Users/robertshear/Projects/netseq-pipeline/test/options.json \
    ~/Projects/netseq-pipeline/netseq-fastq-to-bam.wdl 
