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

        # Genome source for STAR
        File refFasta
        Int MultimapNmax = 1

        # Unprocessed reads
        File inputFastQ
        String sampleName = basename(basename(inputFastQ, ".gz"), ".fastq")

        # environment
        String netseq_docker = 'rdshear/netseq'
        Int preemptible = 1
    }

    call StarAlign {
        input:
            Infile = inputFastQ,
            refFasta = refFasta,
            sampleName = sampleName,
            MultimapNmax = MultimapNmax,
            threads = threads,
            docker = netseq_docker,
            preemptible = preemptible
    }

    call BamToBedgraph {
        input:
            AlignedBamFile = StarAlign.output_bam,
            sampleName = sampleName,
            threads = threads,
            docker = netseq_docker,
            preemptible = preemptible
    }

    output {

        Array[File] starLogs = StarAlign.star_logs
        File output_bam = StarAlign.output_bam

        File dedup_bam = BamToBedgraph.BamFileDeduped
        Array[File] bedgraphs = BamToBedgraph.CoverageBedgraphs
        Array[File] dedup_logs = BamToBedgraph.DedupLogs
    }
}

# TODO rename StarAlign to AlignReads
task StarAlign {
    input {
        File Infile
        File refFasta
        String sampleName
        Int MultimapNmax

        Int threads = 8
        String docker
        Int preemptible
    }

    String bamResultName = "~{sampleName}.aligned.bam"

    command <<<
        set -e

        cp ~{refFasta} ./sacCer3.fa
        samtools faidx sacCer3.fa
        gatk CreateSequenceDictionary -R ./sacCer3.fa

        STAR \
        --runMode genomeGenerate \
        --runThreadN ~{threads} \
        --genomeDir star_work \
        --genomeFastaFiles ./sacCer3.fa \
        --genomeSAindexNbases 10
 
        # force the temp directory to the docker's disks
        tempStarDir=$(mktemp -d)
        # star wants to create the directory itself
        rmdir "$tempStarDir"

        samtools import -0 ~{Infile} | \
        python3 /scripts/ExtractUmi.py /dev/stdin /dev/stdout 6 RX \
        | STAR --runMode alignReads \
            --genomeDir star_work \
            --runThreadN ~{threads} \
            --readFilesIn /dev/stdin \
            --readFilesCommand samtools view \
            --readFilesType SAM SE \
            --outTmpDir "$tempStarDir" \
            --outStd SAM \
            --outFileNamePrefix aligned/~{sampleName}. \
            --outReadsUnmapped None \
            --outFilterMultimapNmax ~{MultimapNmax} \
            --clip3pAdapterSeq ATCTCGTATGCCGTCTTCTGCTTG \
            --clip3pNbases 0 \
            --clip5pNbases 6  \
            --alignIntronMax 1 \
        | samtools sort >  ~{bamResultName}
>>>

    # TODO glob the *.out and/or log files
    output {
        File output_bam = bamResultName
        Array[File] star_logs = glob('aligned/*.out')
    }

    runtime {
        docker: docker
        memory: "8G"
        cpu: threads
        disks: "local-disk 25 SSD"
        preemptible: preemptible
    }
}

task BamToBedgraph {
    input {
        File AlignedBamFile
        String sampleName
        String umi_tag = 'RX'

        Int threads
        String docker
        Int preemptible
    }

    String bamDedupName = "~{sampleName}.dedup.bam"

    command <<<
        df
        set -e

        samtools index ~{AlignedBamFile}
        mamba run -n umi_tools umi_tools dedup -I ~{AlignedBamFile} \
            --output-stats=~{sampleName}.dedup.stats.log \
            --log=~{sampleName}.dedup.log \
            -S ~{bamDedupName} \
            --method unique \
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
        preemptible: preemptible
    }

}
