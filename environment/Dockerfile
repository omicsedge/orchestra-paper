FROM registry.codeocean.com/codeocean/py-r:python3.10.12-R4.2.3-IRkernel-ubuntu22.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get -y upgrade && apt-get install -y \
    cmake \
    apt-utils \
    libz-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    parallel \
    bcftools \
    samtools \
    vcftools \
    && apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN wget https://github.com/MesserLab/SLiM/releases/download/v3.7.1/SLiM.zip && unzip SLiM.zip && rm SLiM.zip
RUN cmake SLiM && make -j"$(nproc)" && install slim eidos /usr/bin && rm -rf SLiM

COPY requirements.txt /
RUN pip3 install -r /requirements.txt
