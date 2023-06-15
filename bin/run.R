#!/usr/bin/env Rscript
library(parallel)
library(flowCore)
library(flowWorkspace)
library(openCyto)

# Wrapper function for the primary entrypoint
source("bin/gate.R")
source("bin/impute.R")
source("bin/process.R")
source("bin/save.R")
source("bin/utils.R")

gsl <- process_study()

catf("Saving gating sets")
save_gating_sets(gsl, ".")
catf("Saving gating sets - Done")

print("DONE")