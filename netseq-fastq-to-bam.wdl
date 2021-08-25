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

    output {
        File adapter_metrix = StarAlign.markAdapterMetrics
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
        Int clippingMinimumLength = 24 # TODO Parameterize

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

time python3 <(cat <<CODE
import pysam
# input: bam file with 3' adapter code location at XT tag
# output: bam file hard clipped with RX tag added, reads without adapters removed
# TODO: parameterize umi_length
umi_length = 6
infile = pysam.AlignmentFile("~{sampleName}.withXTtag.bam", mode="rb", check_sq=False)
outfile = pysam.AlignmentFile("/dev/stdout", mode="w", template=infile)

for r in infile.fetch(until_eof=True):
    xt_tag = [i for i in filter(lambda x: x[0] == 'XT', r.tags)]
    if len(xt_tag) > 0:
        s = r.seq
        q = r.qual
        xt_pos = xt_tag[0][1]
        umi = s[0:umi_length]
        umi_qual = q[0:umi_length]
        # TODO? Remove XT tag bcause the reads are trimmed
        # RX: UMI (possibly corrected), QX: quality score for RX
        # OX: original UMI, BZ quality for original UMI
        r.seq = s[umi_length:xt_pos]
        r.qual = q[umi_length:xt_pos]
        r.tags = r.tags + [('RX', umi)]
        outfile.write(r)

outfile.close()
infile.close()
CODE
 )  \
        | STAR --runMode alignReads \
            --genomeDir /home/star_work \
            --runThreadN 3 \
            --readFilesIn  /dev/stdin \
            --readFilesType SAM SE \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix aligned/~{sampleName}. \
            --limitIObufferSize 200000000 \
            --outReadsUnmapped Fastx \
            --outSAMmultNmax 1 \
            --outSAMattributes All \
            --alignSJoverhangMin 1000 \
            --outFilterMismatchNmax 99

        # TODO Should be tmp?
        rm ~{sampleName}.withXTtag.bam  
        mv aligned/~{sampleName}.Aligned.sortedByCoord.out.bam ~{bamResultName}
    >>>
    # TODO merge with uBAM?
    # TODO glob the *.out and/or log files
    output {
        File markAdapterMetrics = metricsFileName
        File output_bam = bamResultName
    }

    runtime {
        docker: docker
        memory: "8G"
        cpu: 3
    }
}

