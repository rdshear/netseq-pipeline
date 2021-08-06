cd ~/temp
Infile='/Users/robertshear/Projects/netseq-pipeline/test/input/xwt-1.fastq.gz'
tmp=$(mktemp).sam
ubamFileName='wt-1.unaligned.sam'

        picard FastqToSam --FASTQ "${Infile}" \
            --OUTPUT "${tmp}" \
            --SM CPAWT1 \
            --PLATFORM illumina \
            --SORT_ORDER queryname
        #TODO add -M sample.txt which is clipped base read count
		picard MarkIlluminaAdapters -I "${tmp}" \
             -O "${ubamFileName}" \
             --THREE_PRIME_ADAPTER NNNNNNATCTCGTATGCCGTCTTCTGCTTG \
            --ADAPTERS SINGLE_END --FIVE_PRIME_ADAPTER GATCGGAAGAGCACACGTCTGAACTCCAGTCAC

        rm "${tmp}"

