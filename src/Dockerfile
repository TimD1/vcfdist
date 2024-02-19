FROM ubuntu:20.04

# set environment variables
ENV LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    PATH=/usr/local/bin:$PATH \
    DEBIAN_FRONTEND=noninteractive \
    LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH

# install packages
RUN apt-get update --fix-missing && \
    yes | apt-get upgrade && \
    apt-get install -y \
    git \
    make \
    g++ \
    curl \
    wget \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev

# set up HTSlib
WORKDIR /tmp
RUN wget https://github.com/samtools/htslib/releases/download/1.17/htslib-1.17.tar.bz2 && \
    tar -xvf htslib-1.17.tar.bz2
WORKDIR ./htslib-1.17
RUN ./configure --prefix=/usr/local && \
    make && \
    make install

# install bcftools
WORKDIR /tmp
RUN wget https://github.com/samtools/bcftools/releases/download/1.17/bcftools-1.17.tar.bz2 && \
    tar -xvf bcftools-1.17.tar.bz2
WORKDIR ./bcftools-1.17
RUN ./configure --prefix=/usr/local && \
    make && \
    make install

# clone repo
WORKDIR /tmp
RUN git clone https://github.com/TimD1/vcfdist
WORKDIR ./vcfdist/src
RUN make && \
    make install
