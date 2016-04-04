FROM google/cloud-sdk
MAINTAINER Zhuyi Xue <zxue.bcgsc@gmail.com>

# linuxbrew setup part is adopted from
# https://raw.githubusercontent.com/sjackman/docker-linuxbrew/master/linuxbrew-core/Dockerfile,
# which was originally written for ubuntu, but google/cloud-sdk is debian (jessie).

RUN apt-get update && apt-get install locales sudo

RUN localedef -i en_US -f UTF-8 en_US.UTF-8

RUN apt-get install -y curl g++ gawk m4 make patch ruby tcl

RUN useradd -m -s /bin/bash linuxbrew
RUN echo 'linuxbrew ALL=(ALL) NOPASSWD:ALL' >>/etc/sudoers

USER linuxbrew
WORKDIR /home/linuxbrew
ENV HOME /home/linuxbrew
ENV PATH /home/linuxbrew/.linuxbrew/bin:/home/linuxbrew/.linuxbrew/sbin:$PATH
ENV SHELL /bin/bash
RUN yes |ruby -e "$(curl -fsSL https://raw.github.com/Homebrew/linuxbrew/go/install)"
RUN brew doctor || true

RUN sudo apt-get install build-essential
