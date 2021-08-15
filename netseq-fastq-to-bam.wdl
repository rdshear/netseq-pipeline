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

        # general environment
        Int threads = 8
        String genome

        # Genome source for STAR
        File? starReferencesIn
        File refFasta
        File refFastaIndex
        File annotationsGTF
        Int starJunctionReadLength = 50 # maximum number of bases to concatanate between donor and acceptor sides of splice junction 

        # Unprocessed reads
        File inputFastQ
        String sampleName = basename(basename(inputFastQ, ".gz"), ".fastq")

        String netseq_docker = 'rdshear/netseq'
    }
    if (!defined(starReferencesIn)) {
        call StarGenerateReferences { 
            input:
                threads = threads,
                genome = genome,
                ref_fasta = refFasta,
                ref_fasta_index = refFastaIndex,
                annotations_gtf = annotationsGTF,
                read_length = starJunctionReadLength,
                docker = netseq_docker
        }
    }
    File starReferences = select_first([StarGenerateReferences.star_genome_refs_zipped,starReferencesIn,""])

    call fastqToSam {
        input:
            sampleName = sampleName,
            Infile = inputFastQ,
            docker = netseq_docker
    }

    call StarAlign {
        input:
            star_genome_refs_zipped = starReferences,
            refFasta = refFasta,
            ubamFile = fastqToSam.ubamFile,
            sampleName = sampleName,
            threads = threads,
            docker = netseq_docker
    }

    output {
        File reads_unmapped_bams = fastqToSam.ubamFile
        File adapter_metrix = fastqToSam.markAdapterMetrics
        File? starReferencesOut = StarGenerateReferences.star_genome_refs_zipped
        Array[File]? starLogs = StarGenerateReferences.star_logs
        File output_bam = StarAlign.output_bam
    }
}

task StarGenerateReferences {
    input {
        String genome
        File ref_fasta
        File ref_fasta_index
        File annotations_gtf
        Int read_length
        Int threads = 8
        String docker
    }
    String starRefsName = "star-~{genome}-refs.tar.gz"
    command <<<
        set -e

        source activate gatk
    
        mkdir star_work

        STAR \
        --runMode genomeGenerate \
        --runThreadN ~{threads} \
        --genomeDir star_work \
        --genomeFastaFiles ~{ref_fasta} \
        --sjdbGTFfile ~{annotations_gtf} \
        --genomeSAindexNbases 10 \
        --sjdbOverhang ~{read_length} \

        tar -zcvf ~{starRefsName} star_work
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
    >>>

    output {
            File ubamFile = ubamFileName
            File markAdapterMetrics = metricsFileName
    }
    runtime {
        docker: docker
    }
}
task StarAlign {
    input {
        File star_genome_refs_zipped
        File ubamFile
        File refFasta
        String sampleName
        Int clippingMinimumLength = 24 # TODO Parameterize

        Int threads = 8
        String docker
    }

    String bamResultName = "~{sampleName}.aligned.bam"

    command <<<
        set -e

        tar -xvzf ~{star_genome_refs_zipped}

        # TODO: should be "along side" of local refFasta
        # TODO Refactor this line
        if [[ ! -f $(dirname ~{refFasta})/$(basename ~{refFasta} '.fa').dict ]]
        then
            gatk CreateSequenceDictionary -R ~{refFasta}
        fi
        
        # TODO use template to avoid empty "touch" file
        fastqFile=$(mktemp).fastq

        # Drop lines without adapters, as indicated by XT filename or STAR will break?
        #TODO Add UMI extraction
        # TODO 
        cat > temp.R <<CODE
            block_size <- 1000
            infile <- file("stdin")
            open(infile)
            while (length(x <- readLines(infile, n = block_size)) > 0 ) {
                writeLines(x[startsWith(x, "@") | grepl("XT:i:", x, fixed = TRUE)])
            } 
        CODE


        samtools view ~{ubamFile} \
        | Rscript --vanilla temp.R \
        | gatk SamToFastq -I all_xt.sam -F $fastqFile --CLIPPING_ACTION X \
            --CLIPPING_ATTRIBUTE XT \
            --CLIPPING_MIN_LENGTH ~{clippingMinimumLength}

        # TODO add parameter --runthreadN ~{threads} 

        STAR --runMode alignReads \
            --genomeDir star_work \
            --readFilesIn $fastqFile \
            --outSAMtype BAM Unsorted \
            --outFileNamePrefix aligned/~{sampleName}. \
            --outReadsUnmapped None \
            --outSAMattributes All \
            --alignIntronMin 11 \
            --alignIntronMax 5000 \
            --outFilterType BySJout \
            --alignSJoverhangMin 8 \
            --alignSJDBoverhangMin 1 \
            --outFilterMismatchNmax 99 \
            --outSAMattrIHstart 0

        samtools sort -n -O sam  aligned/~{sampleName}.Aligned.out.bam > sorted.sam

        gatk MergeBamAlignment --REFERENCE_SEQUENCE ~{refFasta} \
            --ALIGNED sorted.sam \
            --UNMAPPED_BAM ~{ubamFile} \
            --OUTPUT ~{bamResultName}

    >>>

    # TODO glob the *.out and/or log files
    output {
        File output_bam = bamResultName
    }

    runtime {
        docker: docker
    }
}

