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

        
###################### revert to bwa ... but with gatk
bwa index ../inputs/-1031706207/sacCer3.fa



gatk SamToFastq \
I=xwt-1.aligned.bam \
FASTQ=/dev/stdout \
CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true | \
bwa mem -M -t 7 -p ../inputs/-1031706207/sacCer3.fa /dev/stdin > xwt-1.bwa.bam




 | \
SamToFastq -I xwt-1.aligned.bam -FASTQ /dev/stdout -CLIPPING_ATTRIBUTE XT -CLIPPING_ACTION 2 -INTERLEAVE true -NON_PF true | \
gatk MergeBamAlignment \
ALIGNED_BAM=/dev/stdin \
UNMAPPED_BAM=xwt-1.revertsam.bam \
OUTPUT=xwt-1.piped.bam \
R=../inputs/-1031706207/sacCer3.fa CREATE_INDEX=true ADD_MATE_CIGAR=true \
CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS \
TMP_DIR=/path/shlee