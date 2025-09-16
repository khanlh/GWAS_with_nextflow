FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    wget \
    unzip \
    bcftools \
    r-base \
    tzdata \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages('qqman', repos='https://cloud.r-project.org/')"
# Install PLINK 2.0
WORKDIR /opt
RUN wget https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_latest.zip && \
    unzip plink2_linux_x86_64_latest.zip -d /usr/local/bin && \
    rm plink2_linux_x86_64_latest.zip

