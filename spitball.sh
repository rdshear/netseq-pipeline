samtools view ../input/xwt-1.withUMItag.bam \
| STAR --runMode alignReads \
            --genomeDir ../star_work \
            --runThreadN 3 \
            --readFilesIn  /dev/stdin \
            --readFilesType SAM SE \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix aligned/xwt-1. \
            --outReadsUnmapped Fastx \
            --outSAMmultNmax 1 \
            --outSAMattributes All \
            --alignSJoverhangMin 1000 \
            --outFilterMismatchNmax 99 

