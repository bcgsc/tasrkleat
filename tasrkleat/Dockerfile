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

RUN curl -OL "https://github.com/bcgsc/abyss/releases/download/1.5.2/abyss-1.5.2.tar.gz" \
    && tar zxf abyss-1.5.2.tar.gz \
    && cd abyss-1.5.2 \
    && ./configure --with-mpi=/usr/lib/openmpi \
    && make -j 8 \
    && make install \
    && cd .. && rm -rfv abyss-1.5.2*


# install transabyss
RUN apt-get install -y libxml2-dev libpng12-dev
RUN pip install python-igraph pysam biopython
ENV MACHTYPE=x86_64
ENV PATH=${HOME}/bin/${MACHTYPE}:${PATH}
RUN curl -OL http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat \
    && mv -v blat /usr/local/bin
RUN curl -OL https://github.com/bcgsc/transabyss/archive/1.5.2.zip \
    && unzip 1.5.2
ENV PATH=/transabyss-1.5.2:${PATH}


# install samtools
RUN apt-get install -y ncurses-dev
RUN curl -OL https://sourceforge.net/projects/samtools/files/samtools/0.1.18/samtools-0.1.18.tar.bz2 \
    && tar jxf samtools-0.1.18.tar.bz2 \
    && cd samtools-0.1.18 \
    && make
ENV PATH=/samtools-0.1.18:${PATH}


# install bwa
RUN curl -OL https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.12.tar.bz2 \
    && tar jxf bwa-0.7.12.tar.bz2 \
    && cd bwa-0.7.12 \
    && make
ENV PATH=/bwa-0.7.12:${PATH}

    
# install gmap-gsnap
ENV PATH=${HOME}/local/gmap/bin:${PATH}
RUN wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2014-12-28.tar.gz  \
    && tar zxf gmap-gsnap-2014-12-28.tar.gz \
    && cd gmap-2014-12-28 \
    && ./configure --prefix=${HOME}/local/gmap \
    && make \
    && make install \
    && cd .. && rm -rf gmap-gsnap-2014-12-28.tar.gz gmap-2014-12-28


ENV PROJECT_ID=tasrkleat
ENV _PROJECT_DIR ${HOME}/local/${PROJECT_ID}
ENV PATH=${_PROJECT_DIR}:${PATH}
RUN mkdir -p ${_PROJECT_DIR}
ADD app/*.py ${_PROJECT_DIR}/

CMD ["app.py", "--help"]
