FROM ubuntu:22.04
MAINTAINER jmonlong@ucsc.edu

# Prevent dpkg from trying to ask any questions, ever
ENV DEBIAN_FRONTEND noninteractive
ENV DEBCONF_NONINTERACTIVE_SEEN true

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    wget \
    git \
    gcc \
    make \
    bzip2 \
    tabix \
    r-base r-base-dev \
    libxml2-dev \
    libssl-dev \
    libmariadb-dev \
    libcurl4-openssl-dev \
    apt-transport-https \
    software-properties-common \
    dirmngr \
    gpg-agent \
    libncurses5-dev \
    libncursesw5-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    && rm -rf /var/lib/apt/lists/*

## R with sveval package
ADD install.R /home

RUN R -f /home/install.R

## bcftools
RUN wget --no-check-certificate https://github.com/samtools/bcftools/releases/download/1.17/bcftools-1.17.tar.bz2 && \
        tar -xjf bcftools-1.17.tar.bz2 && \
        cd bcftools-1.17 && \
        ./configure && make && make install && \
        cd .. && rm -rf bcftools-1.17.tar.bz2

## Install pandoc
RUN apt-get update \
    && wget https://github.com/jgm/pandoc/releases/download/1.18/pandoc-1.18-1-amd64.deb \
    && dpkg -i pandoc-1.18-1-amd64.deb \
    && apt-get install -f \
    && rm pandoc-1.18-1-amd64.deb

WORKDIR /opt/scripts
ADD annotate_svs.R /opt/scripts
ADD annotate_denovo_svs.R /opt/scripts
ADD annotate_sv_overlap.R /opt/scripts

WORKDIR /home
