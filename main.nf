#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process immune_flow_gating {
    input:
    path "fcs/"
    path "bin/"

    output:
    path "*"

"""#!/usr/bin/env Rscript
source("bin/run.R")
"""

}

workflow {

    if(params.fcs == false){
        log.info("Must specify 'fcs' parameter")
        exit 1
    }

    // Make a channel with all of the R files from bin/
    Channel
        .fromPath(
            "${workflow.projectDir}/bin/*.R",
            glob: true,
            type: 'file')
        .toSortedList()
        .set { rscripts }

    // Make a channel with all of the FCS files in the input
    Channel
        .fromPath(
            "${params.fcs}".split(",").toList()
        )
        .toSortedList()
        .set { fcs }

    immune_flow_gating(
        fcs,
        rscripts
    )

}