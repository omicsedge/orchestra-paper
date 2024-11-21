## hash:sha256:763a54ec0735d2cd9c6a5bd0ceda6b0db1851e5bbcacc4555572d4fd8c9ce9b3
#FROM registry.codeocean.com/codeocean/ubuntu-desktop:ubuntu-22.04-guacamole-1.5.2

# hash:sha256:65ad6591c08d4907cc2a654dc11b7d8e0f6a04978edf4b02ff27526dfb30c456
FROM registry.codeocean.com/codeocean/py-r:python3.10.12-R4.2.3-IRkernel-ubuntu22.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get -y upgrade && apt-get install -y \
    cmake apt-utils libz-dev libbz2-dev liblzma-dev \
    libcurl4-openssl-dev parallel vcftools && \
    apt-get clean && apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN curl http://benhaller.com/slim/SLiM.zip > SLiM.zip && unzip SLiM.zip && rm SLiM.zip
RUN mkdir build && cd build
RUN cmake ../SLiM && make -j"$(nproc)"

RUN git clone --recurse-submodules https://github.com/samtools/htslib.git
RUN git clone https://github.com/samtools/bcftools.git --branch 1.19
RUN cd bcftools && make && make install

COPY requirements.txt /
RUN pip3 install -r /requirements.txt

COPY postInstall /
RUN /postInstall
