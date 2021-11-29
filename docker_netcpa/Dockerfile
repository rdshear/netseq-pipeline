# /bin/bash -c docker build . -t rdshear/seqcpa
#############################################################
# Changepoint analysis of NET-sdeq data sets
#############################################################

# TODO Add user and working directory

FROM continuumio/miniconda3:latest
LABEL description="Docker image with components for NET-seq pipeline"

# Install mamba for faster installation in the subsequent step
RUN conda install -c conda-forge mamba -y

# setup channels so that bioconductor packages load correctly
RUN conda config --add channels bioconda
RUN conda config --add channels r
RUN conda config --add channels conda-forge
RUN conda config --set channel_priority strict

# load the packages we need
RUN mamba install -c conda-forge r-base=4.1.1 r-devtools r-biocmanager -y
RUN mamba install r-tidyverse -y
RUN mamba install bioconductor-GenomicRanges bioconductor-rtracklayer \
      bioconductor-plyranges -y
RUN mamba install h2o -y

# HACK ... need systemfonts for kabelExtra
#       but kabblExtra source install fails (sep-2021)
RUN mamba install r-systemfonts -y
# Install R packages that are possibly not available via conda
COPY bin/install.R /
RUN mamba run -n base Rscript /install.R

# execution time scripts
ADD scripts/ /scripts/


# Record details of conda installation
RUN mamba env export -n base > base_exported.yml

# install OS patches
RUN  apt update &&  apt upgrade -y

# Copy additional scripts from bin and add to PATH
RUN mkdir /opt/bin
COPY bin/* /opt/bin/
RUN chmod +x /opt/bin/*
ENV PATH="$PATH:/opt/bin/"

CMD ["/bin/bash"]
