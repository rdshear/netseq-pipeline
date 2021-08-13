# /bin/bash -c docker build . -t rdshear/netseq
#############################################################
# Net-seq allignment
#############################################################
FROM broadinstitute/gatk
LABEL Maintainer='Robert Shear <rshear@gmail.com>'
USER root
WORKDIR /gatk
RUN conda install -n gatk -c bioconda star
CMD ["/bin/bash"]
