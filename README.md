# TASRCloud

TASRCloud means Targeted ASsembly on the Cloud.

### Download the tasrcloud image

An image has already been built and named **zyxue/tasrcloud**  on
[Docker Hub](https://hub.docker.com/r/zyxue/tasrcloud/). To download it

    $ docker pull zyxue/tasrcloud:latest


### Launch a container interactively for tasrcloud

    $ docker run --rm -t -i zyxue/tasrcloud /bin/bash

BioBloomTools and ABYSS will be available

    root@<container-id>:/# biobloomcategorizer --help

    root@<container-id>:/# abyss-pe --help

### To build the image yourself

    $ git clone git@github.com:bcgsc/tasrcloud.git

    $ cd tasrcloud

    $ docker build <username>/tasrcloud .

If you see a successfully built message, that means the image has been
successfully built. Check with

    $ docker images

For more information about Docker, please try the
[tutorial](https://docs.docker.com/linux/) and read the
[documentation](https://docs.docker.com/).
