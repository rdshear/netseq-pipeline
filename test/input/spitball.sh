cd ~/temp
sampleName="xwt-1"
Infile='/Users/robertshear/Projects/netseq-pipeline/test/input/xwt-1.fastq.gz'
MinimumReadLength='24'
tmp=$(mktemp).sam
ubamFileName='xwt-1.unaligned.sam'
genomeRef='/n/groups/churchman/rds19/starRefFiles/genome/'

conda activate cpa

picard FastqToSam --FASTQ "${Infile}" \
    --OUTPUT "${tmp}" \
    --SM CPAWT1 \
    --PLATFORM illumina \
    --SORT_ORDER queryname

picard MarkIlluminaAdapters -I "${tmp}" \
        -O "${ubamFileName}" \
        -M mia_metrix.txt \
        --THREE_PRIME_ADAPTER NNNNNNATCTCGTATGCCGTCTTCTGCTTG \
    --ADAPTERS SINGLE_END --FIVE_PRIME_ADAPTER GATCGGAAGAGCACACGTCTGAACTCCAGTCAC

picard SamToFastq -I ${ubamFileName} -F tmp.fastq --CLIPPING_ACTION X \
    --CLIPPING_ATTRIBUTE XT \
    --CLIPPING_MIN_LENGTH ${MinimumReadLength}

conda activate star

STAR --genomeDir ${genomeRef}  --readFilesIn tmp.fastq \
         --outFileNamePrefix aligned/${sampleName}. \
         --outReadsUnmapped Fastx \
        --outSAMattributes All \
         --alignIntronMin 11 \
         --alignIntronMax 5000 \
         --outFilterType BySJout \
         --outFilterMultimapNmax 1 \
         --alignSJoverhangMin 8 \
         --alignSJDBoverhangMin 1 \
         --outFilterMismatchNmax 999 \
         --alignMatesGapMax 2000 \
         --outSAMattrIHstart 0 \
         --outSAMtype BAM SortedByCoordinate