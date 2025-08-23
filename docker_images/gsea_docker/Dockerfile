FROM debian:bookworm-slim

SHELL ["/bin/bash", "-c"]

RUN apt-get update && \
    apt-get install -y \
    build-essential \
    unzip \
    python3-dev \
    python3-pip \
    python3-venv \
    default-jre

# Install all software under /opt/software:
RUN mkdir -p /opt/software && \
    mkdir -p /opt/software/gsea && \
    mkdir -p /opt/software/scripts

# Add the various python requirements and script files
ADD requirements.txt /opt/software/
ADD prep_files.py /opt/software/scripts/
ADD move_final_files.py /opt/software/scripts/
ADD create_rnk_file.py /opt/software/scripts/
ADD GSEA_4.3.2.zip /opt/software/gsea/

RUN unzip /opt/software/gsea/GSEA_4.3.2.zip -d /opt/software/gsea

# Install python deps
RUN python3 -m venv /py_venv && \
  /py_venv/bin/pip install --upgrade pip && \
  /py_venv/bin/pip install -r /opt/software/requirements.txt

ENV PATH="/py_venv/bin:$PATH"