FROM google/cloud-sdk
MAINTAINER Zhuyi Xue <zxue.bcgsc@gmail.com>

# crcmod problem is fixed with the following setup.
RUN apt-get update && apt-get -y install \
    gcc \
    python-dev \
    python-setuptools \
    && easy_install -U pip \
    && pip install -U crcmod colorlog ruffus

RUN apt-get -yf install \
    fastqc

# qualimap bamqc needs bam to be sorted, which adds complexity, ignored for now
# 2016-02-23

# RUN wget "https://bitbucket.org/kokonech/qualimap/downloads/qualimap_v2.2.zip" \
#     && unzip qualimap_v2.2.zip \
#     && cd .. && rm -rfv qualimap_v2.2.zip

# ENV PATH=/qualimap_v2.2:${PATH}

ENV PROJECT_ID=bamqc
ENV PATH=/${PROJECT_ID}:${PATH}
RUN mkdir /${PROJECT_ID}
ADD app/*.py /${PROJECT_ID}/


# CMD ["app.py", \
#      "--input-bam", "path/to/some.bam", \
#     ]

