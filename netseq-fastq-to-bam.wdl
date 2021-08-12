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
        String sampleName = basename(inputFastQ, ".fastq")

#        File refFasta
#        File refFastaIndex
#	    File refDict

        String rnaSeq_docker = 'rdshear/netseq'
    }

    call fastqToSam {
        input:
            sampleName = sampleName,
            Infile = inputFastQ,
            docker = rnaSeq_docker
    }
    
    output {
        File reads_unmapped_bams = fastqToSam.outFile
    }
}

task fastqToSam {
    input {
        String sampleName
        File Infile
        String docker
    }

    String ubamFileName = '~{sampleName}.unaligned.bam'

    command <<<
        gatk FastqToSam --FASTQ ~{Infile} \
            --OUTPUT /dev/stdout \
            --SM  ${sampleName} \
            --PLATFORM illumina \
            --SORT_ORDER queryname \
        | gatk MarkIlluminaAdapters -I /dev/stdin \
                -O ~{ubamFileName} \
                -M mia_metrix.txt \
                --THREE_PRIME_ADAPTER NNNNNNATCTCGTATGCCGTCTTCTGCTTG \
            --ADAPTERS SINGLE_END --FIVE_PRIME_ADAPTER GATCGGAAGAGCACACGTCTGAACTCCAGTCAC
            >>>

    output {
        File outFile = ubamFileName
    }

    runtime {
        docker: docker
    }
}


