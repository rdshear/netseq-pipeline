# /bin/bash -c docker build . -t rdshear/netseq
#############################################################
# Net-seq allignment
#############################################################

# Note using the condaforge/mambaforge image due to trouble getting the conda install to find compatibilities in the bioconda packages
FROM condaforge/mambaforge
LABEL Maintainer='Robert Shear <rshear@gmail.com>'

RUN mamba install -c bioconda samtools gatk4 pysam star

CMD ["/bin/bash"]
