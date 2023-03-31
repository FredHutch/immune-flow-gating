#!/usr/bin/env Rscript

# Wrapper function for the primary entrypoint

# Read in the code defined in the workflow repository
# (note that this is actually being staged to a temporary
# directory in the working folder of a Nextflow task)

source("bin/data.R")
source("bin/custom.R")
source("bin/fetch.R")
source("bin/gate.R")
source("bin/impute.R")
source("bin/process.R")
source("bin/save.R")
source("bin/utils.R")

# Process the study files
process_study(
    study = "STUDY",
    input_dir = "fcs"
)

print("DONE")