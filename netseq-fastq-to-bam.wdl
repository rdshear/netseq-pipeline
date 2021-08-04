version 1.0

workflow NETseq {
        meta {
        description: "NET-seq determine RNAP II occupancy. Pre-processes short reads with UMIs in FASTQ format, removes 3' adapter, aligns with reference genome, and produces bam file as well as BEDgraph file of occupancy at each NT position."
        author: "Robert D. Shear"
        email:  "rshear@gmail.com"
    }

    parameter_meta {
        sample_name: {
            description: "Base name of output files.",
            category: "common"
        }
        reads_unmapped_bams: {
            description: "Unaligned reads in BAM format",
            patterns: ["*.bam"]
        }
    }

    input {
        File fastqFile
#        String genomeDir = '/n/groups/churchman/rds19/starRefFiles/genome/'
        String sampleName
        String rnaSeq_docker = 'rdshear/netseq'
    }

    call fastqToSam {
        input:
            sampleName = sampleName,
            Infile = fastqFile,
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

    String outFileName = '~{sampleName}.unaligned.bam'

    command <<<

        source activate gatk
    
        picard FastqToSam --FASTQ ~{Infile} \
            --OUTPUT ~{outFileName} \
            --SM CPAWT1 \
            --PLATFORM illumina
    >>>

    output {
        File outFile = outFileName
    }

    runtime {
        docker: docker
    }
}


