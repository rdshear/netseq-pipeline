#! /bin/bash
cd /work || return
sampleName="xwt-1"
Infile='/ref/GSE159603/xwt-1.fastq.gz'
MinimumReadLength='24'
#tmp=$(mktemp).sam
ubamFileName='xwt-1.unaligned.sam'
genomeRef='/ref/rds19/starRefFiles/genome/'

gatk FastqToSam --FASTQ "${Infile}" \
    --OUTPUT /dev/stdout \
    --SM  ${sampleName} \
    --PLATFORM illumina \
    --SORT_ORDER queryname \
| gatk MarkIlluminaAdapters -I /dev/stdin \
        -O "${ubamFileName}" \
        -M mia_metrix.txt \
        --THREE_PRIME_ADAPTER NNNNNNATCTCGTATGCCGTCTTCTGCTTG \
    --ADAPTERS SINGLE_END --FIVE_PRIME_ADAPTER GATCGGAAGAGCACACGTCTGAACTCCAGTCAC


gatk SamToFastq -I ${ubamFileName} -F tmp.fastq --CLIPPING_ACTION X \
    --CLIPPING_ATTRIBUTE XT \
    --CLIPPING_MIN_LENGTH ${MinimumReadLength}

threads=3

STAR --genomeDir ${genomeRef} --runthreadN ${threads} --readFilesIn tmp.fastq \
         --outFileNamePrefix aligned/${sampleName}. \
         --outReadsUnmapped None \
        --outSAMattributes All \
         --alignIntronMin 11 \
         --alignIntronMax 5000 \
         --outFilterType BySJout \
         --alignSJoverhangMin 8 \
         --alignSJDBoverhangMin 1 \
         --outFilterMismatchNmax 99 \
         --outSAMattrIHstart 0

# TODO SORT by name

gatk MergeBamAlignment -R /ref/rds19/starRefFiles/genome/../sacCer3.fa \
        --ALIGNED aligned/xwt-1.Aligned.sortedByCoord.out.bam \
        --UNMAPPED_BAM xwt-1.unaligned.sam \
        -O merge_alighments.sam

        
