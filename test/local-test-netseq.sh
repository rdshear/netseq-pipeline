cd ~/temp
conda activate cpa
cromwell run ~/Projects/netseq-pipeline/netseq-fastq-to-bam.wdl  -i ~/Projects/netseq-pipeline/test/inputs.json