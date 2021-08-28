conda activate umi_tools
cd /Users/robertshear/Projects/netseq-pipeline/test/wd/ || exit
# -m 12 ... minimum length 6 bases for umi + 8 bases to align = 14
 cutadapt -a ATCTCGTATGCCGTCTTCTGCTTG --discard-untrimmed -m 14  "${1}"  2> cutadapt.log \
| umi_tools extract --bc-pattern=NNNNNN --log=processed.log \
| STAR --runMode alignReads \
    --genomeDir ../star_work \
    --runThreadN 3 \
    --readFilesIn /dev/stdin \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix aligned/ \
    --outReadsUnmapped Fastx \
    --outSAMmultNmax 1 \
    --outSAMattributes All \
    --alignSJoverhangMin 1000 \
    --outFilterMismatchNmax 99 

mv aligned/Aligned.sortedByCoord.out.bam aligned.bam
time samtools index aligned.bam

time umi_tools dedup -I aligned.bam --output-stats=deduplicated -S deduplicated.bam

time Rscript --vanilla ../../scripts/Dedup.R . deduplicated.bam wt-1
