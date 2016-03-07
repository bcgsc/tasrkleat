# reference: https://hub.docker.com/r/genomicpariscentre/star/~/dockerfile/

FROM google/cloud-sdk
MAINTAINER Zhuyi Xue <zxue@bcgsc.ca>

# check wheter crcmod is installed with "gsutil version -l"
# crcmod problem is fixed with the following setup.
RUN apt-get update && apt-get -y install \
    gcc \
    python-dev \
    python-setuptools \
    && easy_install -U pip \
    && pip install -U crcmod colorlog ruffus


RUN apt-get install -y \
    apt-utils \
    build-essential \
    bwa \
    bzip2 \
    gcc-multilib \
    libncurses5-dev  \
    libncursesw5-dev \
    openjdk-7-jdk \
    openjdk-7-jre \
    picard-tools \
    samtools \
    zlib1g-dev


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
       gcc-multilib

RUN apt-get install -y graphviz

ENV PROJECT_ID=gatk
ENV PATH=/${PROJECT_ID}:${PATH}
RUN mkdir /${PROJECT_ID}
ADD app/*.py /${PROJECT_ID}/


# CMD ["app.py", \
#      "--input-gs-bam", "gs://some-data-bucket/test.bam", \
#      "--input-gs-ref-fa", "gs://path/to/hg19.fa" \
#      "--input-gs-ref-dict", "gs://path/to/hg19.dict" \
#      "--upload-gs-bucket", "gs://some-result-bucket", \
#      "--refresh-token", "/refresh-token/refresh-token" \
#     ]
