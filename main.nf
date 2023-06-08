#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process immune_flow_gating {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    container "${params.container}"

    input:
    path "fcs/"
    path "samplesheet.csv"
    path "bin/"

    output:
    path "gs*", type: 'dir'

"""#!/bin/bash
set -e
run.R 2>&1 | tee immune_flow_gating.log
"""

}

workflow {

    if(params.fcs == false){
        log.info("Must specify 'fcs' parameter")
    }

    if(params.samplesheet == false){
        log.info("Must specify 'samplesheet' parameter")
    }

    if(params.outdir == false){
        log.info("Must specify 'outdir' parameter")
    }
    if(params.outdir == false || params.fcs == false){ exit 1 }

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
        file("${params.samplesheet}", checkIfExists: true),
        rscripts
    )

}