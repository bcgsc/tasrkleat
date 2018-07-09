# Tasrkleat, a pipeline for targeted analysis

### What does it do

The pipeline was designed with three types of analysis in mind:

1. *Targeted cleavage site prediction* in candidate genes, which is the main focus
   of this pipeline
1. *Targeted _de novo_ assembly* of candidate genes
1. *Targeted read alignment* for expression quantification of candidate genes,
   i.e. genes of interest

The commonality among the three is that as they are all targeted analysis of
candidate genes instead of all genes available for a given RNA-Seq dataset. In
addition, targeted cleavage sites prediction (Task 1) depends on the results
from the target _de novo_ assembly (Task 3).

Currently, all three analysis will be conducted by default when running the
pipeline. There is still no option implemented to disable any of them yet.

### Citation

Zhuyi Xue, Rene L. Warren, Ewan A. Gibb, Daniel MacMillan, Johnathan Wong, Readman Chiu, S. Austin Hammond, Chen Yang, Ka Ming Nip, Catherine A. Ennis, Abigail Hahn, Sheila Reynolds, Inanc Birol

**Recurrent tumor-specific regulation of alternative polyadenylation of cancer-related genes**

BMC Genomics (in publication)

DOI : 10.1186/s12864-018-4903-7


### Running environment

The pipeline is desigend mainly for the cloud computing environment, the [Google
Cloud Platform](https://cloud.google.com/) (GCP) in particular, and to be used
in the form of a Docker image. However, in principle, nothing prevents it from
being used without Docker in a non-cloud environment.

<!-- Maybe too advanced for general user, for advanced users, they will figure it out anyway -->

<!-- If you intend to use it without docker, make sure you have all the dependencies -->
<!-- installed properly. Please see the included `Dockerfile` for the needed -->
<!-- dependencies. -->

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

Pre-built Docker images are available at the
[dockerhub](https://hub.docker.com/r/zyxue/tasrkleat/tags/). The tags should
match those at the github repo, except for the *v0* tag, which is used for
testing purpose exclusively, and the *latest* tag, which reflects the
automatically built image from the master branch.


### Use docker image

It is recommended to run the pipeline interactively first to get
familiar with its behavior before scaling up the computation.

Fetch an interactive Docker session

```
# You may or may not need sudo depending on your user group setup
sudo docker run -it --rm \
    -v /path/to/reference:/mnt \
    -v /path/to/reads-data/:/data \
    zyxue/tasrkleat:latest \
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

Once you are inside a tasrkleat container as root. The environment looks like

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
    --input-bf /mnt/targets.bf \
    --transabyss-kmer-sizes 32 52 72 \
    --reference-genome /mnt/hg19.fa \
    --reference-genome-gmap-index /mnt/gmapdb \
    --gtf /mnt/ensembl.fixed.sorted.gz
```

* `--input-bf` accepts the pre-built input bloomfilters.
* `--transabyss-kmer-sizes` accept three kmer sizes.
* `--input-tar` could be a tarball of gzipped fastq files, or a gzipped tar of
   uncompressed fastq files, both situation occurs in the TCGA samples. It's
   dealt in the
   [`extract_tarball`](https://github.com/bcgsc/tasrkleat/blob/master/app/app.py#L26)
   function if you need more details. Currently, tasrkleat can only handle
   paired-end data.

After you get familiar with how the pipeline works, you could run it
in batch mode, e.g.

```
sudo docker run --rm \
	-v /path/to/reference:/mnt \
	-v /path/to/reads-data/:/data \
	zyxue/tasrkleat:latest \
    app.py --input-tar /data/data.tar \
           --input-bf /mnt/targets.bf \
           --transabyss-kmer-sizes 32 52 72 \
           --reference-genome /mnt/hg19.fa \
           --reference-genome-gmap-index /mnt/gmapdb \
           --gtf /mnt/ensembl.fixed.sorted.gz
```

The command is mostly the same to that in the interactive mode except for the
parts that enable interaction (e.g. `-it` and `/bin/bash`) are removed. Now it
runs `app.py` directly instead of `/bin/bash` when the container starts.


### Development

1. Version every package installed if possible in the Dockerfile.
2. Push each versioned image explicitly with `docker push zyxue/tasrkleat:<tag>`.
