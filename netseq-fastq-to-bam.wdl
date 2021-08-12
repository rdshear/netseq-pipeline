version 1.0

workflow NETseq {
        meta {
        description: "NET-seq determine RNAP II occupancy. Pre-processes short reads with UMIs in FASTQ format, removes 3' adapter, aligns with reference genome, and produces bam file as well as BEDgraph file of occupancy at each NT position."
        author: "Robert D. Shear"
        email:  "rshear@gmail.com"
    }

    parameter_meta {
        inputFastQ: "Unprocessed reads"
    }
    input {
        File inputFastQ
        String sampleName = basename(basename(inputFastQ, ".gz"), ".fastq")

        String rnaSeq_docker = 'rdshear/netseq'
    }

    call fastqToSam {
        input:
            sampleName = sampleName,
            Infile = inputFastQ,
            docker = rnaSeq_docker
    }
    
    output {
        File reads_unmapped_bams = fastqToSam.ubamFile
        File adapter_metrix = fastqToSam.markAdapterMetrics
    }
}

task fastqToSam {
    input {
        String sampleName
        File Infile
        String docker
    }

    String tempUbamFileName = '~{sampleName}.raw_unaligned.bam'
    String ubamFileName = '~{sampleName}.unaligned.bam'

    String metricsFileName = '~{sampleName}.markAdaptersMetrics.txt'

    command <<<
        gatk FastqToSam --FASTQ ~{Infile} \
            --OUTPUT ~{tempUbamFileName} \
            --SAMPLE_NAME ~{sampleName} \
            --PLATFORM illumina \
            --SORT_ORDER queryname

        gatk MarkIlluminaAdapters --INPUT ~{tempUbamFileName} \
            --OUTPUT ~{ubamFileName} \
            --METRICS ~{metricsFileName} \
            --THREE_PRIME_ADAPTER NNNNNNATCTCGTATGCCGTCTTCTGCTTG \
            --FIVE_PRIME_ADAPTER GATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
            --ADAPTERS SINGLE_END 

        rm ~{tempUbamFileName}
            >>>

    output {
            File ubamFile = ubamFileName
            File markAdapterMetrics = metricsFileName

    }
    runtime {
        docker: docker
    }
}


