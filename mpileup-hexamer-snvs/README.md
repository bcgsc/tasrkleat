## Versions

* pysam: 0.8.2.1
* pybedtools: 0.6.2
* samtools: 0.1.19
* bedtools: v2.24.0
* vcftools: v0.1.12


## Setup development environment

For the best development experience, you need to have Docker installed already

##### Create virtualenv

```
virtualenv venv
. venv/bin/activate
pip install -r requirements.txt
```

##### Build the image

```
docker build -t <username>/<image-name> .
```

##### Start a Docker container

Open another terminal and start a new container by

```
docker run --rm -ti \
	-v ${PWD}/<sample-data-dir>:/<sample-data-dir> \
	-v ${PWD}/app:/app \
	<username>/<image-name> /bin/bash
```

`<sample-data-dir>` is where you have some small data for testing available, so
the above command mounts that directory under the root directory in the docker
container. Note that a docker container has its own isolated filesystem as well
as network, etc. For more info, see
[Docker security](https://docs.docker.com/engine/security/security/). Besides,
the above command also mounts the `app` directory. Brief explanation of the
flags

* `--rm`: remove the container when it exits
* `-ti`: equivalent to `-t -i`, meaning fetch an interactive TTY
* `-v`: mount a host volume to the container's filesystem `/bin/bash`: start a bash shell after launching the container. /bin/bash will be the only process in the container, you could	confirm this with `top` inside the container.

##### Edit and run

At this point, you should have two terminals, namely

* T1: the original terminal
* T2: the one with a container running interactively 

You then do the edits in T1, and run the program in T2. Consider T2 as the
proper running environment you will use in production. Since `${PWD}/app` is
mounted, the changes should become visible immediately in the container with
command like the following

```
app/app.py \
	--input-data <sample-data-dir>/<some-data-file> \
	--some-other-input <sample-data-dir>/<some-other-input-file>
```

In the case of mpileup-hexamer-snvs in particular, it like

```
app/app.py \
	--input-bam experiment/test.bam \
	--input-ref-fa experiment/hg19.fa \
	--input-ref-gff experiment/utr3.75.gff
```
, assuming `<sample-data-dir>` is named `experiment`.

Have fun!
