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
