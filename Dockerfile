# reference: https://hub.docker.com/r/genomicpariscentre/star/~/dockerfile/

FROM us.gcr.io/isb-cgc-03-0006/gatk
MAINTAINER Zhuyi Xue <zxue@bcgsc.ca>

RUN pip install -U pip
RUN pip install cython

RUN wget https://sourceforge.net/projects/samtools/files/samtools/0.1.19/samtools-0.1.19.tar.bz2 \
    && tar jxf samtools-0.1.19.tar.bz2 \
    && cd samtools-0.1.19 \
    && make \
    && find . -executable -exec cp -v '{}' /usr/local/bin ';' \
    && cd .. && rm -rf samtools-0.1.19.tar.bz2 samtools-0.1.19.tar.bz2

RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.24.0/bedtools-2.24.0.tar.gz \
    && tar zxf bedtools-2.24.0.tar.gz \
    && cd bedtools2 \
    && make \
    && make install \
    && cd .. && rm -rf bedtools2 bedtools-2.24.0.tar.gz

RUN pip install \
    pysam==0.8.2.1 \
    pybedtools==0.6.2

ENV PROJECT_ID=mpileup-hexmer-snv
ENV PATH=/${PROJECT_ID}:${PATH}
RUN mkdir /${PROJECT_ID}
ADD app/*.py /${PROJECT_ID}/


# CMD
# app/app.py --input-bam example.sorted.bam --input-ref-fa hg19.fa --input-ref-gff utr3.75.gff
