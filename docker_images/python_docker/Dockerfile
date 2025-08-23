FROM debian:bookworm-slim

SHELL ["/bin/bash", "-c"]

RUN apt-get update && \
    apt-get install -y \
    build-essential \
    git \
    curl \
    python3-dev \
    python3-pip \
    python3-venv

# Install all software under /opt/software:
RUN mkdir -p /opt/software 
RUN mkdir -p /opt/software/resources
RUN mkdir -p /opt/software/scripts

# Download the prior files into this Docker image:
RUN curl https://granddb.s3.amazonaws.com/optPANDA/ppi/ppi_complete.txt -o /opt/software/resources/tissues_ppi.tsv
RUN curl https://granddb.s3.amazonaws.com/optPANDA/motifs/Hugo_motifCellLine.txt -o /opt/software/resources/tissues_motif.tsv
RUN curl ftp://ftp.broadinstitute.org://pub/gsea/msigdb/human/annotations/Human_Ensembl_Gene_ID_MSigDB.v2023.1.Hs.chip -o /opt/software/resources/Human_Ensembl_Gene_ID_MSigDB.v2023.1.Hs.chip

# Add the various python requirements and script files
ADD requirements.txt /opt/software/

# Install python deps
RUN python3 -m venv /py_venv && \
  /py_venv/bin/pip install --upgrade pip && \
  /py_venv/bin/pip install -r /opt/software/requirements.txt

ENV PATH="/py_venv/bin:$PATH"

# Install the NetZoo python library
RUN cd /opt/software && \
    curl -sL https://github.com/netZoo/netZooPy/archive/refs/tags/0.10.9.tar.gz | tar -xz && \
    /py_venv/bin/pip install -e /opt/software/netZooPy-0.10.9/

ADD map_ensg_to_symbol.py /opt/software/scripts/
ADD overlay_ppi.py /opt/software/scripts/
