# Tasrkleat, a pipeline for targeted analysis

### What does it do

The pipeline was designed with multiple usages in mind:

* Targeted read alignment for expression quantification of candidate genes, i.e. genes of interest
* Targeted _de novo_ assembly of candidate genes
* Targeted cleavage sites predictions in candidate genes

The three are related in a way as they are all targeted studies of candidate
genes only instead of all genes available in a given RNA-Seq dataset. Besides,
targeted cleavage sites detection depends on the results from the target _de
novo_ assembly.

However, currently, all three analysis will be conducted by default when running
the pipeline. There is still no option implemented to disable any of them.

### Running environment

The pipeline was desigend mainly for the cloud computing environment, the
[Google Cloud Platform](https://cloud.google.com/) (GCP) in particular, and to be
used in the form of a Docker image. In principle, nothing prevents it from being
used without Docker in a non-cloud environment. However, there are still
GCP-specific commands (e.g.
[`gsutil`](https://cloud.google.com/storage/docs/gsutil)) inside the pipeline,
so it's **NOT** ready for usage outside of GCP yet.

<!-- The user just needs to make sure all the dependencies are installed properly -->
<!-- wherever it's used (See the included `Dockerfile` for dependencies). -->

### Install

The easiest way to install is to build a Docker image with the included
`Dockerfile`, and use that image directly. To build a Docker image, make sure
you have [Docker](https://www.docker.com/) installed, then try

```
git clone git@github.com:bcgsc/tasrkleat.git
cd tasrkleat
make build
```

To see if the image has been built successfully

```
docker images
```
