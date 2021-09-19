version 1.0

workflow netsq_to_changepoint {
        meta {
        description: "Discover changepoints within genes in NETS-seq occupancy"
        author: "Robert D. Shear"
        email:  "rshear@gmail.com"
    }

    # TODO more parameter metadata
    parameter_meta {
        inputFastQ: "Unprocessed reads"
    }
    input {

        String sampleName
        File genelist
        File CoverageBedgraph_Pos
        File CoverageBedgraph_Neg
        Int ShardCount
        Int MaxGenes
        Int GeneTrimLength

        # environment
        Int threads = 8
        String netcpa_docker = 'rdshear/netcpa'
    }

    String workdir = "./workdir/"

    call CreateShards {
        input:
            genelist = genelist,
            ShardCount = ShardCount,
            MaxGenes = MaxGenes,
            workdir = workdir,
            docker = netcpa_docker
    }
}
    task CreateShards {
        input {
        File genelist
        Int ShardCount
        Int MaxGenes
        String workdir

        String docker
        }

        command <<<
            set -e

            Rscript --vanilla /scripts/DiscoverBreakpointsScatter.R \
             ~{genelist} \
             ~{MaxGenes} \
             ~{ShardCount} \
             ~{workdir}
        >>>
    runtime {
        docker: docker
        }

    

    #${sep=", " array_value}
    }