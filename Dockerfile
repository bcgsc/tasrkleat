FROM google/cloud-sdk
MAINTAINER Zhuyi Xue <zxue.bcgsc@gmail.com>

RUN apt-get update

# mpirun needs ssh to run successfully
# bc is from bsdmainutils, which is needed to run abyss
RUN apt-get -yf install \
    automake \
    bedtools \
    bsdmainutils \
    build-essential \
    curl \
    libboost-all-dev \
    libopenmpi-dev \
    libpython-dev \
    libsparsehash-dev \
    libsqlite3-dev \
    openmpi-bin \
    samtools \
    ssh \
    unzip \
    zlib1g-dev

# enable fast gsutil cp large files
RUN apt-get install -yf gcc python-dev python-setuptools \
    && easy_install -U pip \
    && pip install -U crcmod

# RUN curl -OL "http://www.bcgsc.ca/platform/bioinfo/software/biobloomtools/releases/2.0.12/biobloomtools-2.0.12.tar.gz" \
#     && tar zxf biobloomtools-2.0.12.tar.gz \
#     && cd biobloomtools-2.0.12 \
#     && ./configure \
#     && make \
#     && make install \
#     && cd .. && rm -rfv biobloomtools-2.0.12*

# this version can parse @CO header line successfully
RUN curl -OL "https://github.com/zyxue/biobloom/archive/master.zip" \
    && unzip master.zip && cd biobloom-master \
    && ./autogen.sh \
    && ./configure \
    && make -j 8 \
    && make install \
    && cd .. && rm -rfv master.zip biobloom-master

RUN curl -OL "https://github.com/bcgsc/abyss/releases/download/1.9.0/abyss-1.9.0.tar.gz" \
    && tar zxf abyss-1.9.0.tar.gz \
    && cd abyss-1.9.0 \
    && ./configure --with-mpi=/usr/lib/openmpi \
    && make -j 8 \
    && make install \
    && cd .. && rm -rfv abyss-1.9.0*

# This path will differ if the base image is ubtuntu, i.e.
# ENV PATH=/root/anaconda2/bin:${PATH}
ENV PATH=/tasrcloud:/anaconda2/bin:${PATH}

# gsutils only supports python 2.6.x or 2.7.x at the time 2016-01-22
RUN curl -OL  "https://3230d63b5fc54e62148e-c95ac804525aac4b6dba79b00b39d1d3.ssl.cf1.rackcdn.com/Anaconda2-2.4.1-Linux-x86_64.sh" \
    && bash  Anaconda2-2.4.1-Linux-x86_64.sh -b \
    && rm -v Anaconda2-2.4.1-Linux-x86_64.sh

RUN pip install --upgrade \
    colorlog \
    ruffus

# # && pip install --upgrade google-api-python-client \

RUN mkdir /tasrcloud
ADD *.py /tasrcloud/

# CMD value should be json, must use double quotes
CMD ["app.py", \
     "--input-gs-bam", "gs://tasrcloud-test-data/test.bam", \
     "--input-gs-bf", "gs://tasrcloud-bfs/targetUTRcell2009.tar.gz", \
     "--abyss-kmer-size", "15", \
     "--upload-gs-bucket", "gs://tasrcloud-test-results", \
     "--refresh-token", "/refresh-token/refresh-token" \
    ]
