# reference: https://hub.docker.com/r/genomicpariscentre/star/~/dockerfile/

FROM google/cloud-sdk
MAINTAINER Zhuyi Xue <zxue@bcgsc.ca>

# Update the repository sources list
RUN apt-get update && apt-get install -y \
    apt-utils \
    build-essential \
    bwa \
    bzip2 \
    gcc-multilib \
    libncurses5-dev  \
    libncursesw5-dev \
    openjdk-7-jdk \
    openjdk-7-jre \
    picard \
    samtools \
    zlib1g-dev \
    git

# RUN git clone https://github.com/alexdobin/STAR.git /tmp/STAR
# WORKDIR /tmp/STAR
# RUN git checkout -b 2.5.1b
# WORKDIR /tmp/STAR/source

RUN wget https://github.com/alexdobin/STAR/archive/2.5.1b.tar.gz \
    && tar zxf 2.5.1b.tar.gz \
    && cd STAR-2.5.1b \
    && make STAR \
    && mv -v source/STAR /usr/local/bin/ \
    && cd .. && rm -rfv 2.5.1b.tar.gz STAR-2.5.1b

# install GenomeAnalysisTK.jar to /
RUN wget ftp://ftp.bcgsc.ca/public/zxue/GenomeAnalysisTK-3.5.tar.bz2 \
    && tar jxf GenomeAnalysisTK-3.5.tar.bz2 \
    && rm -v GenomeAnalysisTK-3.5.tar.bz2

# Cleanup                                                                                                                                                                                                                                                                                                             RUN rm -rf /tmp/STAR ; apt-get clean ; apt-get remove --yes --purge build-essential gcc-multilib apt-utils zlib1g-dev vim-common git
RUN rm -rf /tmp/STAR \
    && apt-get clean \
    && apt-get remove -y --purge \
       apt-utils \
       build-essential \
       gcc-multilib \
       git
