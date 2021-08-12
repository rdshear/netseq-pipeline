# /bin/bash -c docker build . -t rdshear/netseq
#############################################################
# Net-seq allignment
#############################################################
FROM broadinstitute/gatk
LABEL Maintainer='Robert Shear <rshear@gmail.com>'
USER root
WORKDIR /gatk
RUN conda create -c bioconda -n gatk star
CMD ["/bin/bash"]
