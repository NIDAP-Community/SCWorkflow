FROM nciccbr/ccbr_ubuntu_22.04:v4

# build time variables
ARG BUILD_DATE="000000"
ENV BUILD_DATE=${BUILD_DATE}
ARG BUILD_TAG="000000"
ENV BUILD_TAG=${BUILD_TAG}
ARG REPONAME="000000"
ENV REPONAME=${REPONAME}

ARG R_VERSION=4.1.3
ENV R_VERSION=${R_VERSION}

SHELL ["/bin/bash", "-lc"]

# Install conda and give write permissions to conda folder
RUN echo 'export PATH=/opt2/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh" -O ~/miniforge3.sh && \
    /bin/bash ~/miniforge3.sh -b -p /opt2/conda && \
    rm ~/miniforge3.sh && chmod 777 -R /opt2/conda/
ENV PATH="/opt2/conda/bin:$PATH"

# Pin channels and update
RUN conda config --add channels conda-forge \
 && conda config --add channels bioconda \
 && conda config --set channel_priority strict

# install conda packages
RUN mamba install -y \
    r-base=${R_VERSION} \
    r-anndata=0.7.5.2 \
    r-biocmanager \
    r-callr=3.7.3 \
    bioconductor-celldex=1.4.0 \
    r-colorspace=2.0-3 \
    bioconductor-complexheatmap=2.10.0 \
    r-cowplot=1.1.1 \
    r-data.table=1.15.4 \
    r-dendextend=1.16.0 \
    r-dendsort=0.3.4 \
    r-digest=0.6.37 \
    r-dplyr=1.1.4 \
    bioconductor-edger=3.36.0 \
    r-future=1.34.0 \
    r-future.apply=1.11.2 \
    r-gargle=1.3.0 \
    r-gdata=2.18.0.1 \
    r-ggextra=0.10.1 \
    r-ggplot2=3.3.6 \
    r-ggpubr=0.4.0 \
    r-ggrepel=0.9.5 \
    r-globals=0.16.3 \
    r-gridbase=0.4-7 \
    r-gridextra=2.3 \
    r-gtable=0.3.5 \
    r-harmony \
    bioconductor-hdf5r=1.3.5 \
    r-htmlwidgets=1.6.4 \
    r-httpuv=1.6.15 \
    r-httr=1.4.7 \
    r-jsonlite=1.8.8 \
    r-leiden=0.4.3 \
    bioconductor-limma \
    r-magrittr=2.0.3 \
    r-markdown=1.13 \
    bioconductor-mast=1.20.0 \
    r-pheatmap=1.0.12 \
    r-plotly=4.10.4 \
    r-plyr=1.8.7 \
    r-png=0.1-7 \
    r-progressr=0.14.0 \
    r-purrr=1.0.2 \
    r-quantmod=0.4.20 \
    r-rcolorbrewer=1.1-3 \
    r-reshape2=1.4.4 \
    r-reticulate=1.40.0 \
    r-rlang=1.1.4 \
    r-scales=1.2.1 \
    bioconductor-scdblfinder=1.8.0 \
    r-seurat=4.1.1 \
    bioconductor-singler=1.8.1 \
    r-statmod=1.5.0 \
    r-stringr=1.5.1 \
    r-svglite=2.1.0 \
    r-tibble=3.2.1 \
    r-tidyr=1.2.1 \
    r-tidyverse=1.3.2 \
    r-viridislite=0.4.1 \
    r-xfun=0.47 \
    r-zip=2.3.3 \
    r-knitr=1.48 \
    r-rmarkdown=2.28 \
    r-roxygen2=7.2.3 \
    r-testthat=3.1.6 \
    r-usethis=3.1.0 \
    r-cffr \
    r-covr \
    r-goodpractice \
    r-here=1.0.1 \
    r-lintr \
    r-pkgdown=2.0.7 \
    r-rcmdcheck=1.4.0 \
  && conda clean -afy

# install R package
COPY . /opt2/SCWorkflow
RUN R -e "devtools::install_local('/opt2/SCWorkflow', dependencies = TRUE, upgrade = 'never', repos='http://cran.rstudio.com')"

# add scworkflow exec to the path
# RUN chmod -R +x /opt2/conda/lib/R/library/SCWorkflow/exec
# ENV PATH="$PATH:/opt2/conda/lib/R/library/SCWorkflow/exec"
# RUN scworkflow --help

# copy example script & json to data
COPY ./inst/extdata/example_script.sh /data2/
COPY ./inst/extdata/json_args/ /data2/json_args/

# Save Dockerfile in the docker
COPY Dockerfile /opt2/Dockerfile_${REPONAME}.${BUILD_TAG}
RUN chmod a+r /opt2/Dockerfile_${REPONAME}.${BUILD_TAG}

# cleanup
WORKDIR /data2
RUN apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*