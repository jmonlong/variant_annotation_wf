FROM ubuntu:20.04
MAINTAINER jmonlong@ucsc.edu

# Prevent dpkg from trying to ask any questions, ever
ENV DEBIAN_FRONTEND noninteractive
ENV DEBCONF_NONINTERACTIVE_SEEN true

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    wget \
    gcc \
    git \
    make \
    bzip2 \
    tabix \
    python3 \
    python3-pip \
    libncurses5-dev \
    libncursesw5-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    autoconf \
    build-essential \
    pkg-config \
    apt-transport-https software-properties-common dirmngr gpg-agent \
    && rm -rf /var/lib/apt/lists/*

## samtools
RUN wget --no-check-certificate https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2 && \
    tar -xjvf samtools-1.17.tar.bz2 && \
    cd samtools-1.17 && \
    ./configure && \
    make && \
    make install

## vg
WORKDIR /bin
RUN wget --no-check-certificate https://github.com/vgteam/vg/releases/download/v1.48.0/vg && \
    chmod +x vg
ENV PATH=$PATH:/bin

## minigraph
RUN git clone https://github.com/lh3/minigraph.git && \
    cd minigraph && \
    make
ENV PATH=$PATH:/bin/minigraph

# python packages
RUN pip3 install cyvcf2

# add script
WORKDIR /opt/scripts
ADD validate-svs.py /opt/scripts

WORKDIR /home
