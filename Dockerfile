# /bin/bash -c docker build . -t rdshear/netseq
#############################################################
# Net-seq allignment
#############################################################

# Note using the condaforge/mambaforge image due to trouble getting the conda install to find compatibilities in the bioconda packages
FROM condaforge/mambaforge
LABEL Maintainer='Robert Shear <rshear@gmail.com>'

# TODO drop pysam and bwa 
RUN mamba install -c bioconda samtools cutadapt gatk4 pysam star bwa

CMD ["/bin/bash"]
