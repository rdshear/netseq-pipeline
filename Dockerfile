# /bin/bash -c docker build . -t netseq
#############################################################
# Net-seq allignment
#############################################################
FROM broadinstitute/gatk
LABEL Maintainer='Robert Shear <rshear@gmail.com>'
USER root
WORKDIR /gatk
RUN /bin/bash -c 'source /gatk/gatkenv.rc'
RUN conda config --add channels bioconda
RUN conda install -c bioconda -n gatk picard

# by default /bin/bash is executed
CMD ["/bin/bash"]   
