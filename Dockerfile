FROM python:3.12-slim

# Add build dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    gcc \
    wget \
    zip \
    cmake \
    libgsl-dev \
    autoconf \
    automake \
    libtool \
    bcftools \
    samtools \
    vcftools \
    r-base \
    python3-dev \
    pkg-config \
    libhdf5-dev \
    libz-dev \
    && apt-get clean && apt-get purge && rm -rf /var/lib/apt/lists/*

RUN wget https://github.com/MesserLab/SLiM/releases/download/v3.7.1/SLiM.zip && unzip SLiM.zip && rm SLiM.zip
RUN cmake SLiM && make -j"$(nproc)" && install slim eidos /usr/bin && rm -rf SLiM

COPY environment/requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Set environment variables to handle ARM64 compilation issues
ENV CFLAGS="-O2"
ENV CXXFLAGS="-O2"
ENV NUMCODECS_DISABLE_AVX2=1

RUN pip install --no-cache-dir -r requirements.txt

COPY environment/main.py .

COPY . .

# Set the entry point to the Python script
ENTRYPOINT ["python", "main.py"]