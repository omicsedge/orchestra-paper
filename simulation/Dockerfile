FROM ubuntu:20.04
WORKDIR /app

# Avoid prompts from apt
ENV DEBIAN_FRONTEND=noninteractive

# Update and install necessary packages
RUN apt-get update
RUN apt-get install -y \
    wget \
    python3.8 \
    python3-pip \
    build-essential \
    cmake \
    bcftools \
    samtools \
    vcftools \
    libgsl-dev \
    autoconf \
    automake \
    libtool

# Install R
RUN apt-get install -y r-base

# Install SLIM and tskit
ADD docker_scripts/install_slim.sh .
RUN bash install_slim.sh

# Requirements
ADD requirements.txt ./
RUN pip install -r ./requirements.txt

# Setup python3
RUN ln -s /usr/bin/python3 /usr/bin/python

# Working directory
COPY . .

ENTRYPOINT ["python", "simulate.py"]
