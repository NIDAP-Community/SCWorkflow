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
# Note: Most version pins removed to allow conda to resolve compatible versions with R 4.1.3
# Only R version is pinned per requirement
RUN mamba install -y \
    r-base=${R_VERSION} \
    r-anndata \
    r-biocmanager \
    r-callr \
    bioconductor-celldex \
    r-colorspace \
    bioconductor-complexheatmap \
    r-cowplot \
    r-data.table \
    r-dendextend \
    r-dendsort \
    r-digest \
    r-dplyr \
    bioconductor-edger \
    r-future \
    r-future.apply \
    r-gargle \
    r-gdata \
    r-ggextra \
    r-ggplot2 \
    r-ggpubr \
    r-ggrepel \
    r-globals \
    r-gridbase \
    r-gridextra \
    r-gtable \
    r-harmony \
    r-hdf5r \
    r-htmlwidgets \
    r-httpuv \
    r-httr \
    r-jsonlite \
    r-leiden \
    bioconductor-limma \
    r-magrittr \
    r-markdown \
    bioconductor-mast \
    r-pheatmap \
    r-plotly \
    r-plyr \
    r-png \
    r-progressr \
    r-purrr \
    r-quantmod \
    r-rcolorbrewer \
    r-reshape2 \
    r-reticulate \
    r-rlang \
    r-scales \
    bioconductor-scdblfinder \
    r-seurat=4.1.1 \
    bioconductor-singler \
    r-statmod \
    r-stringr \
    r-svglite \
    r-tibble \
    r-tidyr \
    r-tidyverse \
    r-viridislite \
    r-xfun \
    r-zip \
    r-knitr \
    r-rmarkdown \
    r-roxygen2 \
    r-testthat \
    r-usethis \
    r-cffr \
    r-covr \
    r-goodpractice \
    r-here \
    r-lintr \
    r-pkgdown \
    r-rcmdcheck \
  && conda clean -afy

# install R package
COPY . /opt2/SCWorkflow
RUN R -e "devtools::install_local('/opt2/SCWorkflow', dependencies = TRUE, upgrade = 'never', repos='http://cran.rstudio.com')"

# add scworkflow exec to the path
# RUN chmod -R +x /opt2/conda/lib/R/library/SCWorkflow/exec
# ENV PATH="$PATH:/opt2/conda/lib/R/library/SCWorkflow/exec"
# RUN scworkflow --help

# copy example script & json to data
# COPY ./inst/extdata/example_script.sh /data2/
# COPY ./inst/extdata/json_args/ /data2/json_args/

# Save Dockerfile in the docker
COPY Dockerfile /opt2/Dockerfile_${REPONAME}.${BUILD_TAG}
RUN chmod a+r /opt2/Dockerfile_${REPONAME}.${BUILD_TAG}

# Verify all dependencies from DESCRIPTION are installed
RUN cat > /tmp/check_description_deps.R << 'EOF'
# Parse DESCRIPTION file and check if all dependencies are installed
desc_file <- "/opt2/SCWorkflow/DESCRIPTION"
if (!file.exists(desc_file)) {
  stop("DESCRIPTION file not found at ", desc_file)
}
# Read and parse DESCRIPTION
desc <- read.dcf(desc_file)
# Extract dependencies
extract_packages <- function(str) {
  if (is.na(str) || str == "") return(character(0))
  # Split by comma and clean up whitespace and version specs
  pkgs <- strsplit(str, ",")[[1]]
  pkgs <- trimws(pkgs)
  pkgs <- gsub("\\s*\\(.*\\)$", "", pkgs)  # Remove version specs
  pkgs <- pkgs[pkgs != ""]
  pkgs
}
deps <- unique(c(
  extract_packages(desc[1, "Depends"]),
  extract_packages(desc[1, "Imports"]),
  extract_packages(desc[1, "Suggests"]),
  extract_packages(desc[1, "Config/Needs/dev"])
))
# Remove base R
deps <- deps[!grepl("^R$", deps)]
# Check if each dependency is installed
missing <- deps[!vapply(deps, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing) > 0) {
  stop("The following dependencies are missing: ", paste(missing, collapse = ", "))
} else {
  message("All dependencies are installed.")
}
EOF
RUN R --vanilla --slave --file=/tmp/check_description_deps.R

# cleanup
WORKDIR /data2
RUN apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*