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

### Pre-built docker images

Pre-built Docker images are available at
https://hub.docker.com/r/zyxue/tasrkleat/tags/. Their tags should
match those at the github repo, except for the *latest* tag, which
reflects the built image from the master branch.

### Use docker image

It is recommended to run the pipeline interactively first to get
familiar with its behavior before scaling up the computation.

Fetch an interactive Docker session

```
# You may or may not need sudo depending on your user group setup
sudo docker run \
	-it \
	--rm \
	-v /path/to/reference:/mnt \
	-v /path/to/reads-data/:/data zyxue/tasrkleat:latest \
	/bin/bash
```

`-it` means fetching an interactive pseudo-tty session. For details of
docker run, please see the
[doc](https://docs.docker.com/engine/reference/run/).

`--rm` means to remove the container after it finishes (e.g. you exit
the container). This is optional, but I find it handy. Otherwise, you
will need to cleanup all the finished container manually with `docker
rm`.

`-v` mounts path of local file system to that inside the container so
that the data is accessible by the pipeline. The above command mounts
two paths, one for the references data, and one for the reads data.

`/bin/bash` means to run `bash` when the container first starts so
that you could interact with it.

`reference` should contains all the necessary reference files, a copy
of those used in the manuscript can be found at
http://bcgsc.ca/downloads/tasrkleat-static/on-cloud/.

Now, you are inside a tasrkleat container as root. The environment looks like

```
root@b7aed8a3b50f:/# whoami
root
```

While tasrkleat Docker image is just a binary file with all necessary
software packaged in, a docker container is an running instance of the
image. In analogy to programming, the image is like a class, and the
container is like an instance of that class.

A example command to run the pipeline inside the container

```
app.py \
	--input-tar /data/data.tar \
	--input-bf /mnt/targets.bf -k 32 52 72 \
	--reference-genome /mnt/hg19.fa \
	--reference-genome-gmap-index /mnt/gmapdb \
	--gtf /mnt/ensembl.fixed.sorted.gz \
	--output-gsc-path gs://tasrkleat-benchmark-kleat"
```

`--input-bf` accepts the pre-built input bloomfilters. `-k` accept
three kmer sizes.  Please please the arguments accordingly.

`data.tar` could be a tarball of gzipped paired-end fastq files, or a
gzipped tar of uncompressed fastq files, both situation occurs in the
TCGA samples. It's dealt in the
[`extract_tarball`](https://github.com/bcgsc/tasrkleat/blob/master/app/app.py#L26)
function if you need more details. It only dealt with paired-end reads.

If you indeed want to upload the results to a Google cloud storage
bucket, please make sure it exists and you have the proper
permission. Otherwise, that step will simply fail currently. But the
results would still be on the local file system.


After you get familiar with how the pipeline works, you could run it
in batch mode, e.g.

```
sudo docker run \
	--rm \
	-v /path/to/reference:/mnt \
	-v /path/to/reads-data/:/data \
	zyxue/tasrkleat:latest \
	app.py \
		--input-tar /data/data.tar \
		--input-bf /mnt/targets.bf -k 32 52 72 \
		--reference-genome /mnt/hg19.fa \
		--reference-genome-gmap-index /mnt/gmapdb \
		--gtf /mnt/ensembl.fixed.sorted.gz \
		--output-gsc-path gs://tasrkleat-benchmark-kleat"; done
```

The command is mostly the same to that in the interactive mode except
for the part that enables interaction (i.e. `-it` and `/bin/bash` are
removed). Now it runs `app.py` directly instead of `/bin/bash` when
the container starts.

### Development

1. Version every package installed if possible in the Dockerfile.
2. Run `docker push zyxue/tasrkleat:<tag>` mannually to avoid accidental build for now.
