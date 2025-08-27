FROM debian:bookworm-slim

RUN apt-get update && \
  apt-get install -y build-essential \
    wget \
    curl \
    zlib1g-dev \
    libcurl4-openssl-dev \
    libncurses-dev \
    libbz2-dev \
    liblzma-dev \
    libxml2-dev \
    libopenblas-dev

# Install all software under /opt/software:
RUN mkdir -p /opt/software 
RUN mkdir -p /opt/software/scripts

# Install various R/bioc packages:
ADD environment.yml /opt/software/

ENV MINIFORGE_VERSION=23.3.1-0
ENV CONDA_DIR=/opt/conda
RUN curl -fsSL https://github.com/conda-forge/miniforge/releases/download/$MINIFORGE_VERSION/Miniforge3-$MINIFORGE_VERSION-$(uname)-$(uname -m).sh -o /tmp/miniforge.sh && \
    /bin/bash /tmp/miniforge.sh -b -p $CONDA_DIR && \
    rm /tmp/miniforge.sh

# Activate conda by default and add to PATH
ENV PATH=$CONDA_DIR/bin:$PATH

# Ensure conda and mamba are up-to-date
RUN conda update -y conda && \
    conda install -y mamba -c conda-forge

# Install packages with mamba
ENV ENV_NAME="r-env"
RUN mamba env create -f /opt/software/environment.yml -n "$ENV_NAME"

ENV PATH="/opt/conda/envs/r-env/bin:${PATH}"

ADD deseq2.R /opt/software/scripts/
ADD download_cohort.R /opt/software/scripts/
ADD download_methylation.R /opt/software/scripts/
ADD run_go.R /opt/software/scripts/