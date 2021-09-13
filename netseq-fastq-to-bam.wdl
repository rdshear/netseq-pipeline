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
        Int MultimapNmax = 1

        # Unprocessed reads
        File inputFastQ
        String sampleName = basename(basename(inputFastQ, ".gz"), ".fastq")

        # environment
        String netseq_docker = 'rdshear/netseq'
        Int preemptible = 1
    }
    if (!defined(starReferencesIn)) {
        call StarGenerateReferences { 
            input:
                threads = threads,
                genome = genome,
                ref_fasta = refFasta,
                ref_fasta_index = refFastaIndex,
                docker = netseq_docker,
                preemptible = preemptible
        }
    }
    File starReferences = select_first([StarGenerateReferences.star_genome_refs_zipped,starReferencesIn,""])

    call StarAlign {
        input:
            Infile = inputFastQ,
            star_genome_refs_zipped = starReferences,
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
        Int threads = 8
        String docker
        Int preemptible
    }
    String starRefsName = "star-~{genome}-refs.tar.gz"
    command <<<
        set -e

        mkdir star_work

        # TODO should be optional parameters (14 is default, 10 recommended for sacCer3)
        # or can compute by size of fasta ... round(log(.) / 2 - 1, 0)
        STAR \
        --runMode genomeGenerate \
        --runThreadN ~{threads} \
        --genomeDir star_work \
        --genomeFastaFiles ~{ref_fasta} \
        --genomeSAindexNbases 10

        tar -zcvf ~{starRefsName} star_work
    >>>

    output {
        Array[File] star_logs = glob("*.out")
        File star_genome_refs_zipped = starRefsName
    }

    runtime {
        docker: docker
        cpu: threads
        preemptible: preemptible
    }
}

# TODO rename StarAlign to AlignReads
task StarAlign {
    input {
        File Infile
        File star_genome_refs_zipped
        File refFasta
        String sampleName
        Int MultimapNmax

        Int threads = 8
        String docker
        Int preemptible
    }

    String bamResultName = "~{sampleName}.aligned.bam"

    command <<<
        df
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
 
        # force the temp directory to the docker's disks
        tempStarDir=$(mktemp -d)
        # star wants to create the directory itself
        rmdir "$tempStarDir"

        samtools import -0 ~{Infile} | \
        python3 /scripts/ExtractUmi.py /dev/stdin /dev/stdout 6 RX \
        | STAR --runMode alignReads \
            --genomeDir /home/star_work \
            --runThreadN ~{threads} \
            --readFilesIn /dev/stdin \
            --readFilesCommand samtools view \
            --readFilesType SAM SE \
            --outTmpDir "$tempStarDir" \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix aligned/~{sampleName}. \
            --outReadsUnmapped Fastx \
            --outFilterMultimapNmax ~{MultimapNmax} \
            --clip3pAdapterSeq ATCTCGTATGCCGTCTTCTGCTTG \
            --clip3pNbases 0 \
            --clip5pNbases 6  \
            --alignIntronMax 1


        mv  aligned/~{sampleName}.Aligned.sortedByCoord.out.bam  \
            ~{sampleName}.aligned.bam
>>>

    # TODO glob the *.out and/or log files
    output {
        File output_bam = bamResultName
        # TODO Add the other log files
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
        df
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
