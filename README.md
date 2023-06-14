# immune-flow-gating
Scripts and workflows for automatically gating immune cells from flow cytometry data

## Analyzing Immune Cell Populations

The Nextflow-based workflow implemented in this repository is intended
to be used for the automated gating of immune cell populations from
flow cytometry data.
Starting from a set of FCS files, this workflow will follow the approach
originally outlined in the [HIPCCyto](https://github.com/RGLab/HIPCCyto/)
library to establish the gates needed to identify different immune cell
subsets.

## Analysis Approach

The analysis follows the set of steps which are executed with the
`process_study` command from the [HIPCCyto](https://github.com/RGLab/HIPCCyto/)
library.
However, in order to provide the user with a greater degree of control
over the parameters which are used in that analysis, the analysis scripts
from the [HIPCCyto](https://github.com/RGLab/HIPCCyto/) library have been
reproduced in this repository.

## Development Guide

With the goal of making the minimal set of changes to scripts from
[HIPCCyto](https://github.com/RGLab/HIPCCyto/), our development approach
is to source the functions from those scripts as part of a primary wrapper
function (`bin/run.R`).

Behind the scenes, this is implemented by:

  1. Making a channel with all of the `*.R` files in `$projectDir/bin/`;
  2. Setting up an input for the primary process which places all of those files in `bin/` for the working directory;
  3. The primary entrypoint script (`run.R`) will use `source()` on all dependencies;
  4. Parameters will be mapped into the scripts using the [`env` input type](https://www.nextflow.io/docs/latest/process.html#input-type-env).

## Channel Map

The default mapping of channels for the required markers is:

| alias | channels |
| --- | --- |
| APC-eFluor 780-A | APC-eFluor780-A |
| eFluor 450-A | eFluor450-A |
