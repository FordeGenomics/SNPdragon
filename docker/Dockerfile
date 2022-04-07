FROM ubuntu:20.04
LABEL maintainer="t.cuddihy1@uq.edu.au"
SHELL [ "/bin/bash", "--login", "-c" ]

COPY ./env_* /tmp/

ENV CONDA_DIR /opt/miniconda3

ENV MINICONDA_VERSION 4.9.2
ADD "https://repo.anaconda.com/miniconda/Miniconda3-py38_${MINICONDA_VERSION}-Linux-x86_64.sh" /tmp/miniconda.sh
RUN ls -la /tmp
RUN chmod +x /tmp/miniconda.sh && \
    /tmp/miniconda.sh -b -p $CONDA_DIR && \
    rm /tmp/miniconda.sh

ENV PATH=$CONDA_DIR/bin:$PATH
RUN echo ". $CONDA_DIR/etc/profile.d/conda.sh" >> ~/.profile
RUN conda init bash

RUN conda update --name base --channel defaults conda
RUN conda env create --file /tmp/env_snp_pipeline_general.yml --force
RUN conda env create --file /tmp/env_snp_pipeline_freebayes.yml --force

RUN conda clean --all --yes
RUN rm -f /tmp/env_*

