version 1.0

## TODO: Add licensing (Broad derivative)
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
        Int threads = 8
        String genome
        # Genome source for STAR
        File refFasta
        File refFastaIndex
        Int? readLength
        File? zippedStarReferences
        File annotationsGTF

        # Unprocessed reads
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

    
	if (!defined(zippedStarReferences)) {

		call StarGenerateReferences { 
			input:
                threads = threads,
                genome = genome,
				ref_fasta = refFasta,
				ref_fasta_index = refFastaIndex,
				annotations_gtf = annotationsGTF,
				read_length = readLength,
				docker = rnaSeq_docker
		}
	} 
    
    File starReferences = select_first([zippedStarReferences,StarGenerateReferences.star_genome_refs_zipped,""])

    output {
        File reads_unmapped_bams = fastqToSam.ubamFile
        File adapter_metrix = fastqToSam.markAdapterMetrics
        File starGeneratedReferences = starReferences
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
        set -e

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

task StarGenerateReferences {
    input {
        String genome
        File ref_fasta
        File ref_fasta_index
        File annotations_gtf
        Int read_length = 100
        Int threads = 8
        String docker
    }
    String starRefsName = "star-~{genome}-refs.tar.gz"
    command <<<
        set -e
        mkdir STAR_WORK

        STAR \
        --runMode genomeGenerate \
        --runThreadN ~{threads} \
        --genomeDir STAR_work \
        --genomeFastaFiles ~{ref_fasta} \
        --sjdbGTFfile ~{annotations_gtf} \
        --sjdbOverhang ~{read_length} \
        --runThreadN ~{threads}

        ls STAR_work

        tar -zcvf ~{starRefsName} STAR_work
    >>>

    output {
        Array[File] star_logs = glob("*.out")
        File star_genome_refs_zipped = starRefsName
    }

    runtime {
        docker: docker
        cpu: threads
    }
}

