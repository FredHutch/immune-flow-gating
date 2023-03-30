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
