FROM python:3.12-slim

# Add build dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    gcc \
    bcftools \
    samtools \
    vcftools \
    libgsl-dev \
    autoconf \
    automake \
    libtool \
    && rm -rf /var/lib/apt/lists/*


COPY requirements.txt .
RUN pip install -r requirements.txt

COPY . .

ENTRYPOINT ["python", "inference.py"]
