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

        # environment
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

    call StarAlign {
        input:
            Infile = inputFastQ,
            star_genome_refs_zipped = starReferences,
            refFasta = refFasta,
            sampleName = sampleName,
            threads = threads,
            docker = netseq_docker
    }

    call BamToBedgraph {
        input:
            AlignedBamFile = StarAlign.output_bam,
            sampleName = sampleName,
            threads = threads,
            docker = netseq_docker
    }

    output {
        File adapter_metrix = StarAlign.markAdapterMetrics
        File? starReferencesOut = StarGenerateReferences.star_genome_refs_zipped

        Array[File] starLogs = StarAlign.star_logs
        File output_bam = StarAlign.output_bam

        File dedup_bam = BamToBedgraph.BamFileDeduped
        Array[File] bedgraphs = BamToBedgraph.CoverageBedgraphs
        Array[File] dedup_logs = BamToBedgraph.DedupLogs
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

        # TODO should be optional parameters (14 is default, 10 recommended for sacCer3)
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

task StarAlign {
    input {
        File Infile
        File star_genome_refs_zipped
        File refFasta
        String sampleName

        Int threads = 8
        String docker
    }

    String ubamFileName = '~{sampleName}.unaligned.bam'
    String metricsFileName = '~{sampleName}.markAdaptersMetrics.txt'

    String bamResultName = "~{sampleName}.aligned.bam"

    command <<<
        set -e

        # HACK temporary workaroud while working locally with docker desktop
        # tar is 25x slower on mounted file system vs local disk,
        #so copy from execution directory to /root
        cp ~{star_genome_refs_zipped} /home/star_index.tar.gz
        time tar -xvzf /home/star_index.tar.gz -C /home/


        # TODO: should be "along side" of local refFasta
        # TODO Refactor this line
        # TODO: should be gzip?
        if [[ ! -f $(dirname ~{refFasta})/$(basename ~{refFasta} '.fa').dict ]]
        then
            gatk CreateSequenceDictionary -R ~{refFasta}
        fi
 
        samtools import -0  ~{Infile} | gatk SortSam --INPUT /dev/stdin --OUTPUT ~{ubamFileName} --SORT_ORDER coordinate

        time gatk MarkIlluminaAdapters --INPUT ~{ubamFileName} \
            --OUTPUT ~{sampleName}.withXTtag.bam \
            --METRICS ~{metricsFileName} \
            --THREE_PRIME_ADAPTER ATCTCGTATGCCGTCTTCTGCTTG \
            --FIVE_PRIME_ADAPTER GATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
            --ADAPTERS SINGLE_END

        python3 /scripts/ExtractUmi.py ~{sampleName}.withXTtag.bam ~{sampleName}.withRXtag.bam 6 RX XT
        
        STAR --runMode alignReads \
            --genomeDir /home/star_work \
            --runThreadN ~{threads} \
            --readFilesIn ~{sampleName}.withRXtag.bam \
            --readFilesCommand samtools view \
            --readFilesType SAM SE \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix aligned/~{sampleName}. \
            --outReadsUnmapped Fastx \
            --outFilterType BySJout \
            --outFilterMultimapNmax 1 \
            --alignSJoverhangMin 8 \
            --alignSJDBoverhangMin 1 \
            --outFilterMismatchNmax 999 \
            --alignMatesGapMax 2000 \
            --outSAMattrIHstart 0

        mv aligned/~{sampleName}.Aligned.sortedByCoord.out.bam ~{sampleName}.aligned.bam
    >>>
    # TODO merge with uBAM?
    # TODO glob the *.out and/or log files
    output {
        File markAdapterMetrics = metricsFileName
        File output_bam = bamResultName
        Array[File] star_logs = glob('aligned/*.out')
    }

    runtime {
        docker: docker
        memory: "8G"
        cpu: threads
    }
}

task BamToBedgraph {
    input {
        File AlignedBamFile
        String sampleName
        String umi_tag = 'RX'
        Int threads
        String docker
    }

    String bamDedupName = "~{sampleName}.dedup.bam"

    command <<<
        set -e

        samtools index ~{AlignedBamFile}
        mamba run -n umi_tools umi_tools dedup -I ~{AlignedBamFile} \
            --output-stats=~{sampleName}.dedup.stats.log \
            --log=~{sampleName}.dedup.log \
            -S ~{bamDedupName} \
        --umi-tag=~{umi_tag} --extract-umi-method=tag

        bedtools genomecov -5 -bg -strand - -ibam ~{bamDedupName} | gzip > ~{sampleName}.pos.bedgraph.gz
        bedtools genomecov -5 -bg -strand + -ibam ~{bamDedupName} | gzip > ~{sampleName}.neg.bedgraph.gz
    >>>

        output {
            File BamFileDeduped = bamDedupName
            Array[File] CoverageBedgraphs = ['~{sampleName}.pos.bedgraph.gz', '~{sampleName}.neg.bedgraph.gz']
            Array[File] DedupLogs = glob('*.log')
    }

    runtime {
        docker: docker
        memory: "8G"
        cpu: threads
    }

}
