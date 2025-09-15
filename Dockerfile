FROM nciccbr/ccbr_ubuntu_22.04:v4

# build time variables
ARG BUILD_DATE="000000"
ENV BUILD_DATE=${BUILD_DATE}
ARG BUILD_TAG="000000"
ENV BUILD_TAG=${BUILD_TAG}
ARG REPONAME="000000"
ENV REPONAME=${REPONAME}

ARG R_VERSION=4.3.2
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
RUN mamba install -y -c conda-forge \
    r-base=${R_VERSION} \
    r-devtools r-testthat \
    r-anndata \
    r-callr r-colorspace r-cowplot \
    r-data.table r-dendextend r-dendsort r-digest r-dplyr \
    r-future r-future.apply \
    r-gargle r-gdata r-ggExtra r-ggplot2 r-ggpubr r-ggrepel r-globals r-glue r-gridBase r-gridExtra r-gtable \
    r-harmony r-hdf5r r-htmlwidgets r-httpuv r-httr \
    r-jsonlite \
    r-leiden \
    r-magrittr r-markdown \
    r-pheatmap r-plotly r-plyr r-png r-progressr r-pryr r-purrr \
    r-quantmod \
    r-RColorBrewer r-reshape2 r-reticulate r-rlang \
    r-scales r-Seurat r-statmod r-stringr r-svglite \
    r-tibble r-tidyr r-tidyverse \
    r-viridisLite \
    r-xfun \
    r-zip \
    bioconductor-celldex bioconductor-ComplexHeatmap \
    bioconductor-edger \
    bioconductor-limma \
    bioconductor-MAST \
    bioconductor-scDblFinder bioconductor-SingleR \
    bioconductor-genomicranges \
    bioconductor-summarizedexperiment \
  && conda clean -afy

# install R package
COPY . /opt2/SCWorkflow
RUN R -e "devtools::install_local('/opt2/SCWorkflow', dependencies = TRUE, repos='http://cran.rstudio.com')"

# add scworkflow exec to the path
RUN chmod -R +x /opt2/conda/lib/R/library/SCWorkflow/exec
ENV PATH="$PATH:/opt2/conda/lib/R/library/SCWorkflow/exec"
RUN scworkflow --help

# copy example script & json to data
COPY ./inst/extdata/TestRunjson.sh /data2/
COPY ./inst/extdata/json_args/ /data2/json_args/

# Save Dockerfile in the docker
COPY Dockerfile /opt2/Dockerfile_${REPONAME}.${BUILD_TAG}
RUN chmod a+r /opt2/Dockerfile_${REPONAME}.${BUILD_TAG}

# cleanup
WORKDIR /data2
RUN apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
