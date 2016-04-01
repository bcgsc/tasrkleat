FROM google/cloud-sdk
MAINTAINER Zhuyi Xue <zxue.bcgsc@gmail.com>

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update

RUN apt-get update && apt-get -y install \
    gcc \
    python-dev \
    python-setuptools \
    && easy_install -U pip \
    && pip install -U crcmod colorlog ruffus

# mpirun needs ssh to run successfully
# bc is from bsdmainutils, which is needed to run abyss
# check wheter crcmod is installed with "gsutil version -l"
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
    zlib1g-dev

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


ENV PROJECT_ID=tasrcloud
ENV PATH=/${PROJECT_ID}:${PATH}
RUN mkdir /${PROJECT_ID}
ADD app/*.py /${PROJECT_ID}/


CMD ["app.py", "--help"]

