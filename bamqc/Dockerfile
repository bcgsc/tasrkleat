FROM google/cloud-sdk
MAINTAINER Zhuyi Xue <zxue.bcgsc@gmail.com>

RUN apt-get update

# from gsult help crcmod:
# To compile and install crcmod:

#   sudo apt-get install gcc python-dev python-setuptools
#   sudo easy_install -U pip
#   sudo pip uninstall crcmod
#   sudo pip install -U crcmod

# but:
# Cannot uninstall requirement crcmod, not installed

# mpirun needs ssh to run successfully
# bc is from bsdmainutils, which is needed to run abyss
# check wheter crcmod is installed with "gsutil version -l"
# crcmod problem is fixed with the following setup.
RUN apt-get -yf install \
    fastqc \
    gcc \
    python-dev \
    python-setuptools \
    && easy_install -U pip \
    && pip install -U crcmod colorlog ruffus

# qualimap bamqc needs bam to be sorted, which adds complexity, ignored for now
# 2016-02-23

# RUN wget "https://bitbucket.org/kokonech/qualimap/downloads/qualimap_v2.2.zip" \
#     && unzip qualimap_v2.2.zip \
#     && cd .. && rm -rfv qualimap_v2.2.zip

# ENV PATH=/qualimap_v2.2:${PATH}

ENV PATH=/bamqc:${PATH}
RUN mkdir /bamqc
ADD *.py /bamqc/

# CMD value should be json, must use double quotes
# example CMD
# CMD ["app.py", \
#      "--input-gs-bam", "gs://some-data-bucket/test.bam", \
#      "--upload-gs-bucket", "gs://some-result-bucket", \
#      "--refresh-token", "/refresh-token/refresh-token" \
#     ]

