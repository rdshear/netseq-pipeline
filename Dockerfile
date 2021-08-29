# /bin/bash -c docker build . -t rdshear/netseq
#############################################################
# Net-seq allignment
#############################################################

# Note using the condaforge/mambaforge image due to trouble getting the conda install to find compatibilities in the bioconda packages
FROM condaforge/mambaforge
LABEL Maintainer='Robert Shear <rshear@gmail.com>'

# TODO Add user and working directory
# TODO drop unused packages
RUN mamba install -c bioconda samtools cutadapt gatk4 pysam star bedtools
RUN mamba create -c bioconda -n umi_tools umi_tools
ADD scripts/ /scripts/


CMD ["/bin/bash"]
