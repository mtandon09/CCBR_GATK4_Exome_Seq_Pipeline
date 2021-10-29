# Base image
FROM ubuntu:20.04

MAINTAINER <vishal.koparde@nih.gov, kuhnsa@nih.gov>

# Create Container filesystem specific 
# working directory and opt directories 
RUN mkdir -p /opt2 && mkdir -p /data2
WORKDIR /opt2 

# Set time zone to US east coast 
ENV TZ=America/New_York
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime \
    && echo $TZ > /etc/timezone

# This section installs system packages required for your project
# If you need extra system packages add them here.
# python/3.8.0
RUN apt-get update \
 && apt-get -y upgrade \
 && DEBIAN_FRONTEND=noninteractive apt-get install -y \
      build-essential \
      bzip2 \
      cmake \
      cpanminus \
      curl \
      g++ \
      gcc \
      gfortran \
      git \
      locales \
      make \
      openjdk-8-jdk \
      parallel \
      pigz \
      python3 \
      python3-pip \
      unzip \
      wget

# Make python3 the default interpreter
# and install Python Packages
RUN ln -s /usr/bin/python3.8 /usr/bin/python
RUN pip3 install --upgrade pip \
	  && pip3 install argparse \
	  && pip3 install numpy \
      && pip3 install pysam \
	  && pip3 install scipy

# Perl fix issue
RUN cpanm FindBin Term::ReadLine

# Set the locale
RUN localedef -i en_US -f UTF-8 en_US.UTF-8

# Other misc dependencies
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y \
      figlet \
      libatlas-base-dev \
      libblas-dev \
      libboost-dev \
      libbz2-dev \
      libexpat1-dev \
      libfreetype6-dev \
      libgd-dev \
      libgd-perl \
      libgs-dev \
      libgsl0-dev \
      libgsl-dev \
      libhtml-template-compiled-perl \
      libicu-dev \
      libjudy-dev \
      liblapack-dev \
      liblzma-dev \
      libmysqlclient-dev \
      libncurses-dev \
      libopenmpi-dev \
      libpng-dev \
      librtmp-dev \
      libssl-dev \
      libssl-dev \
      libtool \
      libxml2-dev \
      libxml-libxml-debugging-perl \
      libxml-opml-simplegen-perl \
      libxslt-dev \
      manpages-dev \
      zlib1g \
      zlib1g-dev \
      zlibc

# Common bioinformatics tools
# bwa/0.7.17  bowtie/1.2.3  bowtie2/2.3.5.1 
# bedtools/2.27.1  bedops/2.4.37  samtools/1.10 
# bcftools/1.10.2  vcftools/0.1.16  trimmomatic/0.39
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y \
      bcftools \
      bedops \
      bedtools \
      bowtie \
      bowtie2 \
      bwa \
      samtools \
      trimmomatic \
      vcftools

# Install GATK/4.2.2
RUN wget https://github.com/broadinstitute/gatk/archive/refs/tags/4.2.1.0.zip \
    && unzip /opt2/4.2.1.0.zip \
    && rm /opt2/4.2.1.0.zip
ENV PATH="/opt2/gatk-4.2.1.0/:$PATH"

# Install dependencies needed to add a new repository over HTTPS
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y \
      gnupg \
      dirmngr \
      ca-certificates \
      apt-transport-https \
      software-properties-common

# Install R (4.0) -- (edgeR needs 3.6) -- and R packages
# ggplot2  dplyr  plotly  htmlwidgets  tidyr and a few extras
# For more information, check out: https://cran.r-project.org/bin/linux/ubuntu/
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
    && add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/' \
    && apt-get -y install r-base r-base-core r-recommended r-base-dev \
    && apt-get -y install libcurl4-openssl-dev libssl-dev libboost-dev libxml2-dev
RUN Rscript -e 'install.packages(c("argparse", "knitr", "tidyverse", "dplyr", "plyr", "plotly", "ggplot2", "RColorBrewer", "htmlwidgets"), repos="http://cran.r-project.org")'
RUN Rscript -e 'install.packages(c("shiny", "gridExtra", "flexdashboard", "rmarkdown", "crosstalk", "DT", "reshape2", "circlize", "viridis"), repos="http://cran.r-project.org")'
RUN Rscript -e 'install.packages("BiocManager"); BiocManager::install(c("limma", "edgeR", "ComplexHeatmap"))'

# Add Dockerfile and argparse.bash script
# and export environment variables
ADD Dockerfile /opt2/base_gatk4_wes.dockerfile
COPY argparse.bash /opt2
RUN chmod -R a+rX /opt2 && chmod -R a+x /opt2/argparse.bash
ENV PATH="/opt2/:$PATH"
ENV TRIMMOJAR="/usr/share/java/trimmomatic-0.39.jar"
WORKDIR /data2

# Clean-up step to reduce size
RUN apt-get clean && apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
