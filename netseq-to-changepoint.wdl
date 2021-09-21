version 1.0

workflow netsq_to_changepoint {

    meta {
        description: "Discover changepoints within genes in NETS-seq occupancy"
        author: "Robert D. Shear"
        email:  "rshear@gmail.com"
    }

    input {

        String sampleName
        String Cp_algorithm
        File genelist
        File CoverageBedgraph_Pos
        File CoverageBedgraph_Neg
        Int ShardCount
        Int MaxGenes
        Int GeneTrimLength

        # environment
        Int threads = 8
        Int preemptible = 1
        String docker_netcpa = 'rdshear/netcpa'
    }

    call CreateShards {
        input:
            genelist = genelist,
            ShardCount = ShardCount,
            MaxGenes = MaxGenes,
            docker = docker_netcpa
    }

    scatter (genespec in CreateShards.shard_specs) {
        String Ofile = 'cp_' + basename(genespec)
        call DiscoverBreakpoints {
            input: 
                genelist = genespec,
                CoverageBedgraph_Pos = CoverageBedgraph_Pos,
                CoverageBedgraph_Neg =CoverageBedgraph_Neg,
                Output_Filename = Ofile,
                MaxGenes = MaxGenes,
                GeneTrimLength = GeneTrimLength,

                docker = docker_netcpa,
                threads = threads,
                preemptible = preemptible
        }
    }

    String results_filename = "~{sampleName}_cp_~{Cp_algorithm}.gff3"

    call GatherShards {
        input:
            outFileName = results_filename,
            shardResults = DiscoverBreakpoints.result_file,
            docker = docker_netcpa
    }
    
    output {
        File changepoint_segments = GatherShards.results
    }
}

task CreateShards {
    input {
    File genelist
    Int ShardCount
    Int MaxGenes

    String docker
    }

    command <<<
        set -e

        Rscript --vanilla /scripts/DiscoverBreakpointsScatter.R \
            ~{genelist} \
            ~{MaxGenes} \
            ~{ShardCount}
    >>>

    runtime {
        docker: docker
    }

    output {
        Array[File] shard_specs = glob("shard_*.gff")
    }
}

task DiscoverBreakpoints {
    input {
        File genelist
        File CoverageBedgraph_Pos
        File CoverageBedgraph_Neg
        String Output_Filename
        Int MaxGenes
        Int GeneTrimLength

        String docker
        Int threads
        Int preemptible
    }

    command <<<
        set -e

        Rscript --vanilla /scripts/DiscoverBreakpointsWorker.R \
            ~{genelist} \
            ~{CoverageBedgraph_Pos} \
            ~{CoverageBedgraph_Neg} \
            ~{Output_Filename} \
            ~{MaxGenes} \
            ~{GeneTrimLength}
    >>>

    output {
        File result_file = Output_Filename
    }

    runtime {
        docker: docker
        preemptible: preemptible
        cpu: threads
    }
}

task GatherShards {
    input {
        String outFileName
        Array[File] shardResults
        String docker
    }

    command <<<
        set -e

        Rscript --vanilla /scripts/DiscoverBreakpointsGather.R \
        ~{outFileName} \
        ~{sep=" " shardResults}
    >>>

    output {
        File results = outFileName
    }

    runtime {
        docker: docker
    }
}