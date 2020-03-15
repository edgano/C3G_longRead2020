FROM debian:jessie
LABEL maintainer="Edgar Garriga <edgano@gmail.com>"

## Docker image as base to build what it's needed
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      ed \
      less \
      vim-tiny \
      wget \
      git \
      python \
      build-essential \
      cmake \
      curl \
      libargtable2-0 \
      ca-certificates \
      libpng12-0 \
      python-biopython \
      python-numpy \ 
      ruby \ 
      python-setuptools \
      default-jdk \
      libpng-dev \
      automake \
      libtool \
      r-base

# install StringTie
RUN git clone https://github.com/gpertea/stringtie && \
  cd stringtie && \
  make release && \
  cp stringtie /bin/ 

# install cufflinks
RUN wget http://ccb.jhu.edu/software/cuffraw/tablemaker-2.1.1.Linux_x86_64.tar.gz && \
  tar -xvf tablemaker-2.1.1.Linux_x86_64.tar.gz && \
  cp tablemaker-2.1.1.Linux_x86_64/tablemaker /bin

#
# Install ballgown 
# https://liorpachter.wordpress.com/2015/08/17/a-sleuth-for-rna-seq/
#
#RUN echo 'source("https://bioconductor.org/biocLite.R")' \
#  && echo 'install.packages("BiocManager")' \
#  && echo 'BiocManager::install("ballgown")' 
RUN echo 'install.packages("BiocManager")' \
&& echo 'library(ballgown)'
