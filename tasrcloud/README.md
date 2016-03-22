# TASRCloud

TASRCloud means Targeted ASsembly on the Cloud.

## Docker

### Download the tasrcloud image

An image has already been built and named **zyxue/tasrcloud**  on
[Docker Hub](https://hub.docker.com/r/zyxue/tasrcloud/). To download it

    $ docker pull zyxue/tasrcloud:latest


### Launch a container interactively for tasrcloud

    $ docker run --rm -t -i zyxue/tasrcloud /bin/bash

### Launch a container with mounted volume

    $ docker run --rm -t -i -v ${PWD}/refresh-token:/refresh-token zyxue/tasrcloud

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

## Google Compute Engine

### Create an compute instance of your choice named dev

    $ gcloud compute instances create dev --zone "us-central1-a" --boot-disk-size 100GB --image container-vm --scopes cloud-platform,storage-rw,bigquery

### ssh to and do whatever you want with it

    $ gcloud compute --project "<project-name>" ssh --zone "us-central1-a" "dev"

## Google Container Engine

### Create a cluster with 3 nodes, each has 2 CPUs & 7.5 GB MEM

    # For a list of machine types: https://cloud.google.com/compute/docs/machine-types#standard_machine_types
    $ gcloud container clusters create tasrcloud --num-nodes 3 --machine-type n1-standard-2

See [here](https://cloud.google.com/compute/docs/machine-types) for other machine types.

### List clusters you've created

    $ gcloud container clusters list

### See the detailed description

    $ gcloud container clusters describe tasrcloud

### Passing cluster credentials to kubectl from Kubernetes

    $ gcloud container clusters get-credentials tasrcloud

### Delete the cluster with CAUTION

    $ gcloud container clusters delete tasrcloud

## Kubernetes

### Create a pod/secret/job/service/replication controller/load balancer

    $ kubectl create -f config.{yaml, json}

### Check logs from a pod

    # content within $() just extracts the pod name.
    $ kubectl logs $(kubectl get pods  --output=jsonpath={.items..metadata.name})

### Accessing Kubernetes UI

Find the server IP, username, password first

    $ kubectl config view | grep 'server\|password\|username'

Then go the the server IP address in a browser and follow the instruction. More
information is available [here](https://github.com/kubernetes/kubernetes/blob/v1.0.6/docs/user-guide/ui.md).
