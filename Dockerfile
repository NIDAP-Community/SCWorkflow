FROM rocker/tidyverse:4.3.2

# build time variables
ARG BUILD_DATE="000000"
ENV BUILD_DATE=${BUILD_DATE}
ARG BUILD_TAG="000000"
ENV BUILD_TAG=${BUILD_TAG}
ARG REPONAME="000000"
ENV REPONAME=${REPONAME}

SHELL ["/bin/bash", "-lc"]

RUN chmod 777 /usr/local/lib/R/site-library /usr/local/lib/R/library

# https://github.com/Bioconductor/bioconductor_docker/blob/7335f85420199679432d2a328c3a59b551b6cfd0/bioc_scripts/install_bioc_sysdeps.sh
RUN apt-get update && apt-get upgrade -y && \
  apt-get install -y --no-install-recommends --allow-unauthenticated \
	automake \
	biber \
	byacc \
	cmake \
	coinor-libcgl-dev \
	coinor-libsymphony-dev \
	coinor-libsymphony-doc \
	curl \
	default-jdk \
	default-libmysqlclient-dev \
	fortran77-compiler \
	ggobi \
	graphviz \
	imagemagick \
	jags \
	libapparmor-dev \
	libarchive-dev \
	libarchive-extract-perl \
	libavfilter-dev \
	libboost-dev \
	libbz2-dev \
	libcairo2-dev \
	libcgi-pm-perl \
	libdbd-mysql-perl \
	libdbi-perl \
	libeigen3-dev \
	libfftw3-dev \
	libfile-copy-recursive-perl \
	libfuse-dev \
	libgdal-dev \
	libgeos-dev \
	libgit2-dev \
	libgl1-mesa-dev \
	libglpk-dev \
	libglu1-mesa-dev \
	libgmp3-dev \
	libgsl0-dev \
	libgslcblas0 \
	libgtk2.0-dev \
	libgtkmm-2.4-dev \
	libhdf5-dev \
	libhdf5-serial-dev \
	libhiredis-dev \
	libjpeg-dev \
	libjpeg-turbo8-dev \
	libjpeg8-dev \
	liblapack-dev \
	liblzma-dev \
	libmagick++-dev \
	libmodule-build-perl \
	libmpfr-dev \
	libmysqlclient-dev \
	libncurses-dev \
	libnetcdf-dev \
	libopenbabel-dev \
	libopenmpi-dev \
	libpcre2-dev \
	libperl-dev \
	libpng-dev \
	libpoppler-cpp-dev \
	libpoppler-glib-dev \
	libpq-dev \
	libproj-dev \
	libprotobuf-dev \
	libprotoc-dev \
	librdf0-dev \
	libreadline-dev \
	librtmp-dev \
	libsasl2-dev \
	libsbml5-dev \
	libssl-dev \
	libtiff5-dev \
	libudunits2-dev \
	libv8-dev \
	libxml-simple-perl \
	libxml2-dev \
	libxpm-dev \
	libxt-dev \
	libz-dev \
	libzmq3-dev \
	mono-runtime \
	mpi-default-bin \
	ocl-icd-opencl-dev \
	openmpi-bin \
	openmpi-common \
	openmpi-doc \
	protobuf-compiler \
	python3-pip \
	sqlite3 \
	tabix \
	tcl8.6-dev \
	tk-dev \
	xfonts-100dpi \
	xfonts-75dpi \
 	liblz4-dev \
    automake \
    cmake \
    default-jre \
    g++ \
    gcc \
    gdb \
    gfortran \
    libcurl4-gnutls-dev \
    make \
    pkg-config

# install R package
COPY . /opt2/SCWorkflow
RUN R -e 'remotes::install_version("Matrix", version="1.6.1"); remotes::install_version("Seurat", version="4.3.0", upgrade="never"); remotes::install_version("SeuratObject", version="4.1.3", upgrade="never")' && \
  R -e "remotes::install_local('/opt2/SCWorkflow', dependencies = TRUE, upgrade='never', repos='http://cran.rstudio.com'); library(SCWorkflow)" && \
  R -s -e "readr::write_tsv(tibble::as_tibble(installed.packages()), '/mnt/r-packages.tsv')"

# add scworkflow exec to the path
RUN chmod -R +x /usr/local/lib/R/site-library/SCWorkflow/exec
ENV PATH="$PATH:/usr/local/lib/R/site-library/SCWorkflow/exec"
RUN scworkflow --help

# Save Dockerfile in the docker
COPY Dockerfile /opt2/Dockerfile_${REPONAME}.${BUILD_TAG}
RUN chmod a+r /opt2/Dockerfile_${REPONAME}.${BUILD_TAG}

# cleanup
WORKDIR /data2
RUN apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
