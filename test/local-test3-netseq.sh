#! /bin/bash
cd ~/temp || exit
rm -R ~/temp/*

conda activate cpa

cromwell run -i ~/Projects/netseq-pipeline/test/inputs_3.json -t wdl \
    -o /Users/robertshear/Projects/netseq-pipeline/test/options.json \
    ~/Projects/netseq-pipeline/netseq-fastq-to-bam.wdl 
