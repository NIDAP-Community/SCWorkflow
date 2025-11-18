FROM nciccbr/ccbr_ubuntu_base_20.04:v7

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
COPY .github/environment.yml /data2/
ENV CONDA_ENV=scw
RUN conda update -n base -c conda-forge conda && \
    mamba env create -n ${CONDA_ENV} -f /data2/environment.yml && \
    echo "conda activate ${CONDA_ENV}" > ~/.bashrc && \
    chmod -R a+rx /opt2
ENV PATH="/opt2/conda/envs/${CONDA_ENV}/bin:$PATH"
ENV R_LIBS_USER="/opt2/conda/envs/${CONDA_ENV}/lib/R/library/"
COPY .github/install.R /data2/

# install R package
COPY . /opt2/SCWorkflow
RUN Rscript /opt2/SCWorkflow/.github/install-pak.R /opt2/SCWorkflow/.github/package-versions.txt && \
	R -e "library(SCWorkflow)" && \
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
