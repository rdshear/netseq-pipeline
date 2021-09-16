# /bin/bash -c docker build . -t rdshear/netseq
#############################################################
# Net-seq allignment
# Ref: https://lifebit.gitbook.io/cloudos/tutorials/combining-dependency-management-with-conda-and-docker
#############################################################

# TODO Add user and working directory

FROM continuumio/miniconda3:latest
LABEL description="Base docker image with conda and util libraries"
ARG ENV_NAME="netseq"

# Install mamba for faster installation in the subsequent step
# Install r-base for being able to run the install.R script
RUN conda install -c conda-forge mamba r-base -y

# Install the conda environment
COPY environment.yml /
RUN mamba env create --quiet --name ${ENV_NAME} --file /environment.yml && conda clean -a

# umi_tools is incompatible with other elements in the base image
RUN mamba create -c bioconda -n umi_tools umi_tools

# Install R packages that are possibly not available via conda
COPY bin/install.R /
RUN Rscript /install.R

# execution time scripts
ADD ../scripts/ /scripts/

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/${ENV_NAME}/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN mamba env export --name ${ENV_NAME} > ${ENV_NAME}_exported.yml

# Copy additional scripts from bin and add to PATH
RUN mkdir /opt/bin
COPY bin/* /opt/bin/
RUN chmod +x /opt/bin/*
ENV PATH="$PATH:/opt/bin/"

CMD ["/bin/bash"]
