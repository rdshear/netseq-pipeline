version 1.0

workflow NETseq {
        meta {
        description: "NET-seq determine RNAP II occupancy. Pre-processes short reads with UMIs in FASTQ format, removes 3' adapter, aligns with reference genome, and produces bam file as well as BEDgraph file of occupancy at each NT position."
        author: "Robert D. Shear"
        email:  "rshear@gmail.com"
    }
    # TODO quality filter
    # TODO more parameter metadata
    parameter_meta {
        # STAR index input
        refFasta: "Genome Reference File, FASTA format"
        # STAR alignment parameters
        inputFastQ: "Illumina Read file, FASTQ format"
        sampleName: "Sample name. If not specified, taken as base name of fastq input file"
        outSAMmultNmax: "The number of alignments returned for each read. If 1, then no multimmappers are returned."
        OutFilterMultiMax: "If a read has multimappers in excess of this paramter, then the read is disreagarded. Defaults"

        #Outputs
        output_bam: "TODO"
        dedup_bam: "TODO"
        bedgraph_pos: "TODO"
        bedgraph_neg: "TODO"

        # Environment
        netseq_docker: "TODO"
        threads: "Number of CPUs to request for task runtimes"
        memory: "TODO"
        preemptible: "TODO"
    }
    input {


        # Genome source for STAR
        File refFasta
        Int outSAMmultNmax = 1     # Default to outputting primary alignment only. (Multimap count still available at )
        Int OutFilterMultiMax = 10 # Default to dropping reads with more than 10 alignments

        # Unprocessed reads
        File inputFastQ
        # TODO...optionally pull off the .1 suffix
        String sampleName = basename(basename(inputFastQ, ".gz"), ".fastq")

        # environment
        String netseq_docker = 'rdshear/netseq'
        Int preemptible = 1
        String memory = "8G"
        Int threads = 8
    }



    call StarAlign {
        input:
            Infile = inputFastQ,
            refFasta = refFasta,
            sampleName = sampleName,
            outSAMmultNmax = outSAMmultNmax,
            OutFilterMultiMax = OutFilterMultiMax,
            threads = threads,
            docker = netseq_docker,
            memory = memory,
            preemptible = preemptible
    }

    call BamToBedgraph {
        input:
            AlignedBamFile = StarAlign.output_bam,
            sampleName = sampleName,
            threads = threads,
            docker = netseq_docker,
            memory = memory,
            preemptible = preemptible
    }

    # TODO clean up intermediate files.
    # TODO output bam files should be optional
    output {
        File output_bam = StarAlign.output_bam

        File dedup_bam = BamToBedgraph.BamFileDeduped
        File bedgraph_pos = BamToBedgraph.CoverageBedgraph_Pos
        File bedgraph_neg = BamToBedgraph.CoverageBedgraph_Neg

        Array[File] logs = [StarAlign.star_log,
                                    StarAlign.star_log_std,
                                    StarAlign.star_log_final,
                                    BamToBedgraph.DedupLogs]
    }
}

# TODO rename StarAlign to AlignReads
task StarAlign {
    input {
        File Infile
        File refFasta
        String sampleName
        Int outSAMmultNmax
        Int OutFilterMultiMax

        Int threads = 8
        String docker
        Int preemptible
        String memory
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
            --outFileNamePrefix ~{sampleName}. \
            --outReadsUnmapped None \
            --outSAMmultNmax ~{outSAMmultNmax} \
            --outFilterMultimapNmax ~{OutFilterMultiMax} \
            --clip3pAdapterSeq ATCTCGTATGCCGTCTTCTGCTTG \
            --clip3pNbases 0 \
            --clip5pNbases 6  \
            --alignIntronMax 1 \
        | samtools sort >  ~{bamResultName}
    >>>

    output {
        File output_bam = bamResultName
        File star_log_final = "~{sampleName}.Log.final.out"
        File star_log = "~{sampleName}.Log.out"
        File star_log_std =  "~{sampleName}.Log.std.out"
    }

    runtime {
        docker: docker
        memory: memory
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
        String memory
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

        bedtools genomecov -5 -bg -strand - -ibam ~{bamDedupName} | bgzip > ~{sampleName}.pos.bedgraph.gz
        bedtools genomecov -5 -bg -strand + -ibam ~{bamDedupName} | bgzip > ~{sampleName}.neg.bedgraph.gz
    >>>

        output {
            File BamFileDeduped = bamDedupName
            File CoverageBedgraph_Pos = '~{sampleName}.pos.bedgraph.gz'
            File CoverageBedgraph_Neg = '~{sampleName}.neg.bedgraph.gz'
            File DedupLogs = '~{sampleName}.dedup.log'
    }

    runtime {
        docker: docker
        memory: memory
        cpu: threads
        preemptible: preemptible
    }

}
