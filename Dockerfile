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
    gcc \
    libboost-all-dev \
    libopenmpi-dev \
    libpython-dev \
    libsparsehash-dev \
    libsqlite3-dev \
    openmpi-bin \
    python-dev \
    python-setuptools \
    samtools \
    ssh \
    unzip \
    zlib1g-dev \
    && easy_install -U pip \
    && pip install -U crcmod colorlog ruffus

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

ENV PATH=/tasrcloud:${PATH}

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

# # && pip install --upgrade google-api-python-client \
