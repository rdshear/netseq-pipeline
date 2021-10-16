# /bin/bash -c docker build . -t rdshear/netseq
#############################################################
# Net-seq allignment
#############################################################

# TODO Add user and working directory

FROM continuumio/miniconda3:latest
LABEL description="Docker image with components for NET-seq pipeline"

# Install mamba for faster installation in the subsequent step
RUN conda install -c conda-forge mamba -y

# Install the needed apps
RUN conda install -c bioconda  sra-tools samtools gatk4 pysam star bedtools


# umi_tools is incompatible with other elements in the base image
RUN mamba create -c bioconda -n umi_tools umi_tools && conda clean -a

# execution time scripts
ADD scripts/ /scripts/

# Dump the details of the installed packages to a file for posterity
RUN mamba env export -n base > base_exported.yml
RUN mamba env export -n umi_tools > umi_tools_exported.yml

CMD ["/bin/bash"]
