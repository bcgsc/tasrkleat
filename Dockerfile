FROM google/cloud-sdk
MAINTAINER Zhuyi Xue <zxue.bcgsc@gmail.com>

RUN apt-get update
RUN apt-get -yf install build-essential automake curl unzip zlib1g-dev libpython-dev

# mpirun needs ssh to run successfully
# bc is from bsdmainutils, which is needed to run abyss
RUN apt-get -yf install \
    libboost-all-dev openmpi-bin \
    libsparsehash-dev libopenmpi-dev libsqlite3-dev ssh bsdmainutils

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
    && ./configure \
    && make -j 8 \
    && make install \
    && cd .. && rm -rfv abyss-1.9.0*

# This path will differ if the base image is ubtuntu, i.e.
# ENV PATH=/root/anaconda2/bin:${PATH}
ENV PATH=/anaconda2/bin:${PATH}

# gsutils only supports python 2.6.x or 2.7.x at the time 2016-01-22
RUN curl -OL  "https://3230d63b5fc54e62148e-c95ac804525aac4b6dba79b00b39d1d3.ssl.cf1.rackcdn.com/Anaconda2-2.4.1-Linux-x86_64.sh" \
    && bash  Anaconda2-2.4.1-Linux-x86_64.sh -b \
    && rm -v Anaconda2-2.4.1-Linux-x86_64.sh

RUN conda create -p /root/venv ipython

RUN /bin/bash -c "source activate /root/venv \
                  && pip install --upgrade pip \
                  && pip install ruffus \
		  && pip install colorlog"
# && pip install --upgrade google-api-python-client \

RUN apt-get install -yf bedtools samtools

RUN mkdir /tasrcloud
ADD *.py run_tasrcloud /tasrcloud/
ENV PATH=/tasrcloud:${PATH}

CMD ["run_tasrcloud"]
