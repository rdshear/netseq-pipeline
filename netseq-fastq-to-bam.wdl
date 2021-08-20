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

        tar -xvzf ~{star_genome_refs_zipped}

        # TODO: should be "along side" of local refFasta
        # TODO Refactor this line
        if [[ ! -f $(dirname ~{refFasta})/$(basename ~{refFasta} '.fa').dict ]]
        then
            gatk CreateSequenceDictionary -R ~{refFasta}
        fi
 
        gatk FastqToSam --FASTQ ~{Infile} \
            --OUTPUT /dev/stdout \
            --SAMPLE_NAME ~{sampleName} \
            --PLATFORM illumina \
            --SORT_ORDER queryname \
        | gatk MarkIlluminaAdapters --INPUT /dev/stdin \
            --OUTPUT /dev/stdout \
            --METRICS ~{metricsFileName} \
            --THREE_PRIME_ADAPTER NNNNNNATCTCGTATGCCGTCTTCTGCTTG \
            --FIVE_PRIME_ADAPTER GATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
            --ADAPTERS SINGLE_END \
| python3 <(cat <<CODE
import pysam

umi_length = 6
base_complements = ''.maketrans({'A':'T', 'C':'G', 'G':'C', 'T':'A'})
infile = pysam.AlignmentFile("/dev/stdin", mode="rb", check_sq=False)
outfile = pysam.AlignmentFile("/dev/stdout", mode="wb", template=infile)

for r in infile.fetch(until_eof=True):
    xt_tag = [i for i in filter(lambda x: x[0] == 'XT', r.tags)]
    if len(xt_tag) > 0:
        # STAR translates seq to reverse complement, restore it to get true UMI
        if r.is_reverse:
            s = r.seq[::-1].translate(base_complements)
            q = r.qual[::-1]
        else:
            s = r.seq
            q = r.qual
        xt_pos = xt_tag[0][1]
        umi = s[xt_pos - 1:xt_pos - 1 + umi_length]
        umi_qual = q[xt_pos - 1:xt_pos - 1 + umi_length]
        # RX: UMI (possibly corrected), QX: quality score for RX
        # OX: original UMI, BZ quality for original UMI
        r.tags = r.tags + [('RX', umi)]
    outfile.write(r)

outfile.close()
infile.close()
CODE
) \
        | gatk SamToFastq -I /dev/stdin -F /dev/stdin --CLIPPING_ACTION X \
            --CLIPPING_ATTRIBUTE XT \
           --CLIPPING_MIN_LENGTH ~{clippingMinimumLength}


        STAR --runMode alignReads \
            --genomeDir star_work \
            --runThreadN ~{threads} \
            --readFilesIn $fastqFile \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix aligned/~{sampleName}. \
            --outReadsUnmapped None \
            --outSAMmultNmax 1 \
            --outSAMattributes All \
            --alignIntronMin 11 \
            --alignIntronMax 5000 \
            --outFilterType BySJout \
            --alignSJoverhangMin 8 \
            --alignSJDBoverhangMin 1 \
            --outFilterMismatchNmax 99 \
            --outSAMattrIHstart 0

        mv aligned/~{sampleName}.Aligned.sortedByCoord.out.bam ~{bamResultName}
    >>>

    # TODO glob the *.out and/or log files
    output {
        File markAdapterMetrics = metricsFileName
        File output_bam = bamResultName
    }

    runtime {
        docker: docker
        memory: "8G"
        cpu: 2

    }
}

