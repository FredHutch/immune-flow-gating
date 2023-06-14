#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process immune_flow_gating {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { it.replaceAll(/.*\//, "") }
    container "${params.container}"

    input:
    path "fcs/"
    path "samplesheet.csv"
    path "bin/"
    path "channel_map.json"

    output:
    path "gs*/gs/*.h5", emit: 'h5'
    path "gs*/gs/*.pb", emit: 'pb'
    path "immune_flow_gating.log", emit: 'log'

"""#!/bin/bash
set -e
run.R 2>&1 | tee immune_flow_gating.log
"""

}

process extract_h5 {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    container "${params.container__python}"

    input:
        path h5

    output:
        path "*.csv.gz"

    """#!/usr/bin/env python3
import h5py
import pandas as pd

# Read the file
input_fp = "${h5}"
print(f"Reading {input_fp}")
dat = h5py.File(input_fp)

# Build a DataFrame
print("Building a DataFrame")
df = pd.DataFrame(
    dat['data'],
    index=[i[0].decode() for i in dat['params']]
).T
print(f"Built a DataFrame with {df.shape[0]:,} rows and {df.shape[1]:,} columns")

# Write to CSV
output_fp = input_fp.replace(".h5", ".csv.gz")
print(f"Writing to {output_fp}")
df.to_csv(output_fp, index=None)
print("Done")
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

    // Get the file with the channel mappings
    channel_map = file(params.channel_map, checkIfExists: true)

    // Make a channel with all of the FCS files in the input
    Channel
        .fromPath(
            "${params.fcs}".split(",").toList(),
            checkIfExists: true
        )
        .toSortedList()
        .set { fcs }

    immune_flow_gating(
        fcs,
        file("${params.samplesheet}", checkIfExists: true),
        rscripts,
        channel_map
    )

    extract_h5(
        immune_flow_gating
            .out
            .h5
            .flatten()
    )

}