#! /bin/bash
cd ~/temp || exit
conda activate cpa

cromwell run -i ~/Projects/netseq-pipeline/test/inputs.json ~/Projects/netseq-pipeline/netseq-fastq-to-bam.wdl 

#cromwell run -o /Users/robertshear/Projects/netseq-pipeline/macos.conf \
#    -i ~/Projects/netseq-pipeline/test/inputs.json \
#     ~/Projects/netseq-pipeline/netseq-fastq-to-bam.wdl 

# java -Dconfig.file=/Users/robertshear/Projects/netseq-pipeline/macos.conf -jar /Users/robertshear/mambaforge/envs/cpa/share/cromwell/cromwell.jar run /Users/robertshear/Projects/netseq-pipeline/netseq-fastq-to-bam.wdl --inputs /Users/robertshear/Projects/netseq-pipeline/test/inputs.json
