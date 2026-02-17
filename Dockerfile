FROM rocker/rstudio:4.1.3

RUN apt-get update && apt-get install -y --no-install-recommends \
      bash gcc g++ gfortran make cmake pkg-config \
      git vim-tiny \
      libcurl4-openssl-dev libssl-dev libxml2-dev \
      libpng-tools \
      libpng-dev libjpeg-dev libtiff5-dev zlib1g-dev \
      libfreetype6-dev libharfbuzz-dev libfribidi-dev \
      libbz2-dev liblzma-dev \
      libgsl-dev libhdf5-dev \
      tcl tk tcl-dev tk-dev \
      libglpk-dev libglpk40 \
      libgit2-dev \
      python3 python3-venv \
      python3-dev python3-pip \
      libfontconfig1-dev \
      autoconf automake libtool \
      libgeos-dev libproj-dev \
      libcairo2-dev libxt-dev \
      automake \
    && rm -rf /var/lib/apt/lists/*


#      libpng libtiff-4 libjpeg libwebp libwebpmux 

# (Optional) prevent R networking during local install steps
# xxx rpf get rid of
# RUN echo 'options(repos = c())' >> /usr/local/lib/R/etc/Rprofile.site

#RUN pip3 install --no-cache-dir igraph leidenalg numpy
# weird pip instal bug introduced 
RUN pip3 install --no-cache-dir igraph "leidenalg==0.10.0" numpy

WORKDIR /home/rstudio
# Copy your local tarballs (ensure PACKAGESS is in the build context)

COPY PACKAGESS/ /opt/pkgs/

# Remove any spatstat* that might be preinstalled in the rocker image
RUN R -q -e " \
  for (lib in .libPaths()) { \
    ip <- rownames(installed.packages(lib.loc = lib)); \
    pk <- grep('^spatstat', ip, value = TRUE); \
    if (length(pk)) { \
      message('Removing from ', lib, ': ', paste(pk, collapse=', ')); \
      remove.packages(pk, lib = lib); \
    } \
  }"


RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/bitops_1.0-8.tar.gz
RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/RCurl_1.98-1.14.tar.gz

RUN set -eux; \
    R CMD INSTALL -l /usr/local/lib/R/site-library \
 /opt/pkgs/GlobalOptions_0.1.2.tar.gz \
 /opt/pkgs/bit_4.0.4.tar.gz \
 /opt/pkgs/BiocGenerics_0.40.0.tar.gz 

#RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/RCurl_1.98.1.16.tar.gz
#RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/RCurl_1.98.1.12.tar.gz


RUN set -eux; \
    R CMD INSTALL -l /usr/local/lib/R/site-library \
 /opt/pkgs/assertthat_0.2.1.tar.gz \
 /opt/pkgs/Biobase_2.54.0.tar.gz \
 /opt/pkgs/bit64_4.0.5.tar.gz \
 /opt/pkgs/shape_1.4.6.tar.gz

RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/remotes_2.5.0.tar.gz
RUN R -e "remotes::install_github('NIDAP-Community/SCWorkflow', ref = 'GalaxyCLI')"

RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/colorspace_2.0-3.tar.gz
RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/circlize_0.4.15.tar.gz

RUN set -eux; \
    R CMD INSTALL -l /usr/local/lib/R/site-library \
 /opt/pkgs/rlang_1.1.4.tar.gz \
 /opt/pkgs/cli_3.6.3.tar.gz \
 /opt/pkgs/glue_1.7.0.tar.gz 

RUN set -eux; \
    R CMD INSTALL -l /usr/local/lib/R/site-library \
 /opt/pkgs/lifecycle_1.0.4.tar.gz \
 /opt/pkgs/vctrs_0.6.5.tar.gz \
 /opt/pkgs/blob_1.2.3.tar.gz \
 /opt/pkgs/clue_0.3-61.tar.gz \
 /opt/pkgs/cluster_2.1.2.tar.gz \
 /opt/pkgs/codetools_0.2-18.tar.gz \
 /opt/pkgs/crayon_1.5.3.tar.gz \
 /opt/pkgs/data.table_1.15.4.tar.gz \
 /opt/pkgs/DBI_1.2.3.tar.gz

RUN set -eux; \
    R CMD INSTALL -l /usr/local/lib/R/site-library \
 /opt/pkgs/dendsort_0.3.4.tar.gz \
 /opt/pkgs/digest_0.6.37.tar.gz \
 /opt/pkgs/ellipsis_0.3.2.tar.gz \
 /opt/pkgs/evaluate_0.24.0.tar.gz \
 /opt/pkgs/fansi_1.0.6.tar.gz \
 /opt/pkgs/farver_2.1.1.tar.gz \
 /opt/pkgs/fastmap_1.2.0.tar.gz \
 /opt/pkgs/cachem_1.1.0.tar.gz \
 /opt/pkgs/fastmatch_1.1-3.tar.gz \
 /opt/pkgs/iterators_1.0.14.tar.gz \
 /opt/pkgs/foreach_1.5.2.tar.gz \
 /opt/pkgs/formatR_1.14.tar.gz

RUN set -eux; \
    R CMD INSTALL -l /usr/local/lib/R/site-library \
 /opt/pkgs/generics_0.1.3.tar.gz \
 /opt/pkgs/rjson_0.2.21.tar.gz \
 /opt/pkgs/GetoptLong_1.0.5.tar.gz \
 /opt/pkgs/gtable_0.3.5.tar.gz \
 /opt/pkgs/labeling_0.4.2.tar.gz \
 /opt/pkgs/munsell_0.5.0.tar.gz \
 /opt/pkgs/R6_2.5.1.tar.gz \
 /opt/pkgs/gridExtra_2.3.tar.gz \
 /opt/pkgs/RColorBrewer_1.1-3.tar.gz \
 /opt/pkgs/utf8_1.2.4.tar.gz \
 /opt/pkgs/pillar_1.9.0.tar.gz 


RUN set -eux; \
    R CMD INSTALL -l /usr/local/lib/R/site-library \
 /opt/pkgs/pkgconfig_2.0.3.tar.gz \
 /opt/pkgs/viridisLite_0.4.1.tar.gz \
 /opt/pkgs/scales_1.2.1.tar.gz \
 /opt/pkgs/withr_3.0.1.tar.gz \
 /opt/pkgs/gtools_3.9.5.tar.gz \
 /opt/pkgs/gridGraphics_0.5-1.tar.gz \
 /opt/pkgs/hms_1.1.2.tar.gz \
 /opt/pkgs/S4Vectors_0.32.4.tar.gz \
 /opt/pkgs/IRanges_2.28.0.tar.gz \
 /opt/pkgs/irlba_2.3.5.1.tar.gz \
 /opt/pkgs/jsonlite_1.8.8.tar.gz 

# /opt/pkgs/jsonlite_1.8.8.tar.gz

RUN set -eux; \
    R CMD INSTALL -l /usr/local/lib/R/site-library \
 /opt/pkgs/png_0.1-7.tar.gz \
 /opt/pkgs/KernSmooth_2.23-20.tar.gz  \
 /opt/pkgs/lambda.r_1.2.4.tar.gz \
 /opt/pkgs/lattice_0.20-45.tar.gz  \
 /opt/pkgs/lazyeval_0.2.2.tar.gz \
 /opt/pkgs/limma_3.50.3.tar.gz \
 /opt/pkgs/locfit_1.5-9.9.tar.gz \
 /opt/pkgs/Matrix_1.5-1.tar.gz \
 /opt/pkgs/matrixStats_0.62.0.tar.gz \
 /opt/pkgs/MatrixGenerics_1.6.0.tar.gz \
 /opt/pkgs/memoise_2.0.1.tar.gz \
 /opt/pkgs/mgcv_1.8-39.tar.gz

RUN set -eux; \
    R CMD INSTALL -l /usr/local/lib/R/site-library \
 /opt/pkgs/nlme_3.1-155.tar.gz \
 /opt/pkgs/rsvd_1.0.5.tar.gz \
 /opt/pkgs/DelayedArray_0.20.0.tar.gz \
 /opt/pkgs/Rcpp_1.0.13.tar.gz \
 /opt/pkgs/sparseMatrixStats_1.6.0.tar.gz \
 /opt/pkgs/DelayedMatrixStats_1.16.0.tar.gz \
 /opt/pkgs/ScaledMatrix_1.2.0.tar.gz \
 /opt/pkgs/uuid_1.1-0.tar.gz \
 /opt/pkgs/xfun_0.47.tar.gz \
 /opt/pkgs/xtable_1.8-4.tar.gz \
 /opt/pkgs/zlibbioc_1.40.0.tar.gz \
 /opt/pkgs/XVector_0.34.0.tar.gz

RUN set -eux; \
    R CMD INSTALL -l /usr/local/lib/R/site-library \
 /opt/pkgs/yaml_2.3.10.tar.gz \
 /opt/pkgs/GenomeInfoDbData_1.2.7.tar.gz \
 /opt/pkgs/GenomeInfoDb_1.30.1.tar.gz \
 /opt/pkgs/beachmat_2.10.0.tar.gz \
 /opt/pkgs/doParallel_1.0.17.tar.gz \
 /opt/pkgs/edgeR_3.36.0.tar.gz \
 /opt/pkgs/futile.options_1.0.1.tar.gz \
 /opt/pkgs/futile.logger_1.4.3.tar.gz \
 /opt/pkgs/GenomicRanges_1.46.1.tar.gz \
 /opt/pkgs/base64enc_0.1-3.tar.gz \
 /opt/pkgs/htmltools_0.5.8.1.tar.gz \
 /opt/pkgs/rappdirs_0.3.3.tar.gz

#RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/tinytex_0.45.tar.gz

RUN set -eux; \
    R CMD INSTALL -l /usr/local/lib/R/site-library \
 /opt/pkgs/jquerylib_0.1.4.tar.gz \
 /opt/pkgs/tinytex_0.44.tar.gz \
 /opt/pkgs/fs_1.6.4.tar.gz \
 /opt/pkgs/sass_0.4.6.tar.gz \
 /opt/pkgs/mime_0.12.tar.gz \
 /opt/pkgs/bslib_0.4.2.tar.gz \
 /opt/pkgs/magrittr_2.0.3.tar.gz \
 /opt/pkgs/stringi_1.8.4.tar.gz \
 /opt/pkgs/stringr_1.5.1.tar.gz \
 /opt/pkgs/highr_0.10.tar.gz \
 /opt/pkgs/knitr_1.48.tar.gz \
 /opt/pkgs/fontawesome_0.5.1.tar.gz

RUN set -eux; \
    R CMD INSTALL -l /usr/local/lib/R/site-library \
 /opt/pkgs/rmarkdown_2.28.tar.gz \
 /opt/pkgs/htmlwidgets_1.6.4.tar.gz \
 /opt/pkgs/sys_3.4.1.tar.gz \
 /opt/pkgs/askpass_1.1.tar.gz \
 /opt/pkgs/openssl_2.0.5.tar.gz \
 /opt/pkgs/curl_6.4.0.tar.gz \
 /opt/pkgs/httr_1.4.7.tar.gz \
 /opt/pkgs/Biostrings_2.62.0.tar.gz \
 /opt/pkgs/cpp11_0.4.7.tar.gz \
 /opt/pkgs/timechange_0.2.0.tar.gz


RUN set -eux; \
    R CMD INSTALL -l /usr/local/lib/R/site-library \
 /opt/pkgs/lubridate_1.8.0.tar.gz \
 /opt/pkgs/isoband_0.2.7.tar.gz \
 /opt/pkgs/tibble_3.2.1.tar.gz \
 /opt/pkgs/ggplot2_3.3.6.tar.gz \
 /opt/pkgs/patchwork_1.2.0.tar.gz \
 /opt/pkgs/later_1.3.2.tar.gz \
 /opt/pkgs/promises_1.3.0.tar.gz \
 /opt/pkgs/crosstalk_1.2.0.tar.gz \
 /opt/pkgs/purrr_1.0.2.tar.gz \
 /opt/pkgs/tidyselect_1.2.1.tar.gz \
 /opt/pkgs/dplyr_1.1.4.tar.gz \
 /opt/pkgs/tidyr_1.2.1.tar.gz

RUN set -eux; \
    R CMD INSTALL -l /usr/local/lib/R/site-library \
 /opt/pkgs/plotly_4.10.4.tar.gz \
 /opt/pkgs/plyr_1.8.7.tar.gz \
 /opt/pkgs/clipr_0.8.0.tar.gz \
 /opt/pkgs/prettyunits_1.1.1.tar.gz \
 /opt/pkgs/progress_1.2.2.tar.gz \
 /opt/pkgs/tzdb_0.4.0.tar.gz \
 /opt/pkgs/vroom_1.6.3.tar.gz \
 /opt/pkgs/readr_2.1.2.tar.gz \
 /opt/pkgs/reshape2_1.4.4.tar.gz \
 /opt/pkgs/plogr_0.2.0.tar.gz \
 /opt/pkgs/RSQLite_2.3.9.tar.gz \
 /opt/pkgs/systemfonts_1.2.3.tar.gz


# /opt/pkgs/vroom_1.6.1.tar.gz \

RUN set -eux; \
    R CMD INSTALL -l /usr/local/lib/R/site-library \
 /opt/pkgs/textshaping_0.3.6.tar.gz \
 /opt/pkgs/ragg_1.2.5.tar.gz \
 /opt/pkgs/dbplyr_2.2.1.tar.gz \
 /opt/pkgs/rstudioapi_0.14.tar.gz \
 /opt/pkgs/dtplyr_1.3.1.tar.gz \
 /opt/pkgs/backports_1.4.1.tar.gz \
 /opt/pkgs/broom_1.0.1.tar.gz \
 /opt/pkgs/ps_1.8.1.tar.gz \
 /opt/pkgs/processx_3.8.5.tar.gz \
 /opt/pkgs/callr_3.7.3.tar.gz \
 /opt/pkgs/reprex_2.0.2.tar.gz \
 /opt/pkgs/modelr_0.1.9.tar.gz

RUN set -eux; \
    R CMD INSTALL -l /usr/local/lib/R/site-library \
 /opt/pkgs/conflicted_1.2.0.tar.gz \
 /opt/pkgs/rematch2_2.1.2.tar.gz \
 /opt/pkgs/gargle_1.3.0.tar.gz \
 /opt/pkgs/rematch_1.0.1.tar.gz \
 /opt/pkgs/cellranger_1.1.0.tar.gz \
 /opt/pkgs/ids_1.0.1.tar.gz \
 /opt/pkgs/googledrive_2.0.0.tar.gz \
 /opt/pkgs/googlesheets4_1.0.1.tar.gz \
 /opt/pkgs/readxl_1.4.1.tar.gz \
 /opt/pkgs/selectr_0.4-2.tar.gz \
 /opt/pkgs/xml2_1.3.3.tar.gz \
 /opt/pkgs/rvest_1.0.3.tar.gz \
 /opt/pkgs/forcats_0.5.2.tar.gz \
 /opt/pkgs/haven_2.5.1.tar.gz \
 /opt/pkgs/tidyverse_1.3.2.tar.gz 

RUN set -eux; \
    R CMD INSTALL -l /usr/local/lib/R/site-library \
 /opt/pkgs/KEGGREST_1.34.0.tar.gz \
 /opt/pkgs/AnnotationDbi_1.56.2.tar.gz \
 /opt/pkgs/snow_0.4-4.tar.gz \
 /opt/pkgs/BH_1.81.0-1.tar.gz \
 /opt/pkgs/BiocParallel_1.28.3.tar.gz \
 /opt/pkgs/BiocSingular_1.10.0.tar.gz \
 /opt/pkgs/ComplexHeatmap_2.10.0.tar.gz \
 /opt/pkgs/viridis_0.6.5.tar.gz


RUN set -eux; \
    R CMD INSTALL -l /usr/local/lib/R/site-library \
 /opt/pkgs/dendextend_1.16.0.tar.gz \
 /opt/pkgs/fgsea_1.20.0.tar.gz \
 /opt/pkgs/ggrepel_0.9.5.tar.gz \
 /opt/pkgs/l2p_0.0-13.tar.gz  \
 /opt/pkgs/l2psupp_0.0-13.tar.gz  \
 /opt/pkgs/ica_1.0-3.tar.gz \
 /opt/pkgs/Rtsne_0.16.tar.gz \
 /opt/pkgs/ggridges_0.5.3.tar.gz \
 /opt/pkgs/scattermore_1.2.tar.gz \
 /opt/pkgs/listenv_0.9.1.tar.gz \
 /opt/pkgs/globals_0.16.3.tar.gz \
 /opt/pkgs/parallelly_1.38.0.tar.gz 


RUN set -eux; \
    R CMD INSTALL -l /usr/local/lib/R/site-library \
 /opt/pkgs/future_1.34.0.tar.gz \
 /opt/pkgs/future.apply_1.11.2.tar.gz \
 /opt/pkgs/RcppEigen_0.3.3.9.3.tar.gz \
 /opt/pkgs/RcppAnnoy_0.0.19.tar.gz \
 /opt/pkgs/zoo_1.8-12.tar.gz \
 /opt/pkgs/lmtest_0.9-40.tar.gz \
 /opt/pkgs/fitdistrplus_1.1-8.tar.gz 

 # /opt/pkgs/gplots_3.1.3.tar.gz \

RUN set -eux; \
    R CMD INSTALL -l /usr/local/lib/R/site-library \
 /opt/pkgs/caTools_1.18.2.tar.gz \
 /opt/pkgs/gplots_3.1.3.tar.gz \
 /opt/pkgs/ROCR_1.0-11.tar.gz \
 /opt/pkgs/igraph_2.0.3.tar.gz \
 /opt/pkgs/pbapply_1.5-0.tar.gz \
 /opt/pkgs/commonmark_1.9.0.tar.gz \
 /opt/pkgs/httpuv_1.6.15.tar.gz \
 /opt/pkgs/sourcetools_0.1.7-1.tar.gz \
 /opt/pkgs/shiny_1.9.1.tar.gz \
 /opt/pkgs/miniUI_0.1.1.1.tar.gz

#RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/igraph_2.0.3.tar.gz


RUN set -eux; \
    R CMD INSTALL -l /usr/local/lib/R/site-library \
 /opt/pkgs/progressr_0.14.0.tar.gz \
 /opt/pkgs/sitmo_2.0.2.tar.gz \
 /opt/pkgs/dqrng_0.4.1.tar.gz \
 /opt/pkgs/FNN_1.1.3.2.tar.gz \
 /opt/pkgs/RcppProgress_0.4.2.tar.gz \
 /opt/pkgs/uwot_0.1.14.tar.gz \
 /opt/pkgs/cowplot_1.1.1.tar.gz 

RUN set -eux; \
    R CMD INSTALL -l /usr/local/lib/R/site-library \
 /opt/pkgs/RcppTOML_0.2.2.tar.gz \
 /opt/pkgs/rprojroot_2.0.4.tar.gz \
 /opt/pkgs/here_1.0.1.tar.gz \
 /opt/pkgs/reticulate_1.40.0.tar.gz \
 /opt/pkgs/leiden_0.4.3.tar.gz \
 /opt/pkgs/RANN_2.6.1.tar.gz \
 /opt/pkgs/RcppArmadillo_0.12.4.0.0.tar.gz \
 /opt/pkgs/sctransform_0.3.4.tar.gz

# see https://github.com/satijalab/seurat/issues/8948
#devtools::install_version("spatstat.core", version = "2.4-0", repos='http://cran.us.r-project.org')
#devtools::install_version("spatstat.data", version = "2.1-2", repos='http://cran.us.r-project.org')
#devtools::install_version("spatstat.geom", version = "2.3-2", repos='http://cran.us.r-project.org')
#devtools::install_version("spatstat.random", version = "2.1-0", repos='http://cran.us.r-project.org')

# old - phils spec ...
# /opt/pkgs/spatstat.core_2.4-4.tar.gz \
# /opt/pkgs/spatstat.data_3.1-2.tar.gz \
# /opt/pkgs/spatstat.geom_3.3-2.tar.gz \
# /opt/pkgs/spatstat.random_3.3-1.tar.gz \
# /opt/pkgs/spatstat.explore_3.1-0.tar.gz \
# /opt/pkgs/spatstat.univar_3.0-0.tar.gz \

#      /opt/pkgs/spatstat.geom_3.2-1.tar.gz \
#      /opt/pkgs/spatstat.data_3.1-2.tar.gz \
#      /opt/pkgs/spatstat.data_2.1-2.tar.gz \
#      /opt/pkgs/spatstat.random_3.1-5.tar.gz \
#      /opt/pkgs/spatstat.geom_2.3-2.tar.gz \

RUN set -eux; \
    R CMD INSTALL -l /usr/local/lib/R/site-library \
      /opt/pkgs/spatstat.utils_3.1-0.tar.gz \
      /opt/pkgs/spatstat.data_3.0-0.tar.gz \
      /opt/pkgs/deldir_2.0-4.tar.gz \
      /opt/pkgs/polyclip_1.10-7.tar.gz \
      /opt/pkgs/spatstat.univar_2.0-3.tar.gz \
      /opt/pkgs/spatstat.geom_3.0-3.tar.gz \
      /opt/pkgs/spatstat.random_2.2-0.tar.gz \
      /opt/pkgs/abind_1.4-5.tar.gz \
      /opt/pkgs/tensor_1.5.tar.gz \
      /opt/pkgs/goftest_1.2-3.tar.gz \
      /opt/pkgs/spatstat.sparse_3.1-0.tar.gz 

#      /opt/pkgs/spatstat.explore_3.0-3.tar.gz \

RUN set -eux; \
    R CMD INSTALL -l /usr/local/lib/R/site-library \
      /opt/pkgs/spatstat.core_2.4-2.tar.gz \
      /opt/pkgs/spatstat.linnet_2.2-1.tar.gz 

RUN set -eux; \
    R CMD INSTALL -l /usr/local/lib/R/site-library \
 /opt/pkgs/RSpectra_0.16-1.tar.gz \
 /opt/pkgs/dotCall64_1.1-1.tar.gz \
 /opt/pkgs/spam_2.11-0.tar.gz  \
 /opt/pkgs/RcppHNSW_0.4.1.tar.gz \
 /opt/pkgs/leidenbase_0.1.30.tar.gz 


RUN set -eux; \
    R CMD INSTALL -l /usr/local/lib/R/site-library \
 /opt/pkgs/fastDummies_1.7.3.tar.gz \
 /opt/pkgs/sp_1.5-0.tar.gz \
 /opt/pkgs/rgeos_0.5-9.tar.gz \
 /opt/pkgs/SeuratObject_4.1.1.tar.gz \
 /opt/pkgs/Seurat_4.1.1.tar.gz \
 /opt/pkgs/XML_3.99-0.14.tar.gz \
 /opt/pkgs/anndata_0.7.5.2.tar.gz \
 /opt/pkgs/annotate_1.72.0.tar.gz \
 /opt/pkgs/ash_1.0-15.tar.gz 


RUN set -eux; \
    R CMD INSTALL -l /usr/local/lib/R/site-library \
 /opt/pkgs/beeswarm_0.4.0.tar.gz \
 /opt/pkgs/filelock_1.0.3.tar.gz \
 /opt/pkgs/BiocFileCache_2.2.1.tar.gz \
 /opt/pkgs/BiocManager_1.30.26.tar.gz \
 /opt/pkgs/BiocNeighbors_1.12.0.tar.gz \
 /opt/pkgs/BiocVersion_3.14.0.tar.gz \
 /opt/pkgs/bluster_1.4.0.tar.gz \
 /opt/pkgs/bslib_0.5.0.tar.gz \
 /opt/pkgs/cachem_1.1.0.tar.gz \
 /opt/pkgs/Cairo_1.6-0.tar.gz \
 /opt/pkgs/desc_1.4.2.tar.gz \
 /opt/pkgs/pkgload_1.3.0.tar.gz \
 /opt/pkgs/brio_1.1.3.tar.gz

#RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/BiocManager_1.30.26.tar.gz

RUN set -eux; \
    R CMD INSTALL -l /usr/local/lib/R/site-library \
 /opt/pkgs/zip_2.3.3.tar.gz \
 /opt/pkgs/openxlsx_4.2.5.2.tar.gz \
 /opt/pkgs/praise_1.0.0.tar.gz \
 /opt/pkgs/diffobj_0.3.5.tar.gz \
 /opt/pkgs/waldo_0.5.1.tar.gz \
 /opt/pkgs/testthat_3.1.6.tar.gz \
 /opt/pkgs/nloptr_2.0.3.tar.gz \
 /opt/pkgs/minqa_1.2.5.tar.gz \
 /opt/pkgs/lme4_1.1-33.tar.gz \
 /opt/pkgs/MatrixModels_0.5-1.tar.gz \
 /opt/pkgs/SparseM_1.81.tar.gz \
 /opt/pkgs/quantreg_5.95.tar.gz \
 /opt/pkgs/numDeriv_2016.8-1.1.tar.gz \
 /opt/pkgs/pbkrtest_0.5.2.tar.gz \
 /opt/pkgs/maptools_1.1-7.tar.gz \
 /opt/pkgs/carData_3.0-5.tar.gz \
 /opt/pkgs/rio_0.5.29.tar.gz \
 /opt/pkgs/car_3.0-0.tar.gz 


RUN set -eux; \
    R CMD INSTALL -l /usr/local/lib/R/site-library \
 /opt/pkgs/pROC_1.18.2.tar.gz \
 /opt/pkgs/proxy_0.4-27.tar.gz \
 /opt/pkgs/e1071_1.7-13.tar.gz \
 /opt/pkgs/ModelMetrics_1.2.2.2.tar.gz \
 /opt/pkgs/SummarizedExperiment_1.24.0.tar.gz \
 /opt/pkgs/SingleCellExperiment_1.16.0.tar.gz 

RUN set -eux; \
    R CMD INSTALL -l /usr/local/lib/R/site-library \
 /opt/pkgs/clock_0.7.0.tar.gz \
 /opt/pkgs/shinyjs_2.1.0.tar.gz \
 /opt/pkgs/colourpicker_1.2.0.tar.gz \
 /opt/pkgs/combinat_0.0-8.tar.gz \
 /opt/pkgs/conquer_1.3.3.tar.gz \
 /opt/pkgs/corrplot_0.92.tar.gz \
 /opt/pkgs/cpp11_0.4.7.tar.gz \
 /opt/pkgs/diagram_1.6.5.tar.gz \
 /opt/pkgs/DT_0.28.tar.gz \
 /opt/pkgs/Rttf2pt1_1.3.12.tar.gz \
 /opt/pkgs/extrafontdb_1.0.tar.gz \
 /opt/pkgs/extrafont_0.18.tar.gz \
 /opt/pkgs/vipor_0.4.7.tar.gz \
 /opt/pkgs/ggbeeswarm_0.7.2.tar.gz \
 /opt/pkgs/ggrastr_1.0.2.tar.gz \
 /opt/pkgs/fastICA_1.2-3.tar.gz \
 /opt/pkgs/gdata_2.18.0.1.tar.gz \
 /opt/pkgs/ggExtra_0.10.1.tar.gz

RUN set -eux; \
    R CMD INSTALL -l /usr/local/lib/R/site-library \
 /opt/pkgs/ggsci_3.0.0.tar.gz \
 /opt/pkgs/ggsignif_0.6.3.tar.gz \
 /opt/pkgs/gower_1.0.1.tar.gz \
 /opt/pkgs/graph_1.72.0.tar.gz \
 /opt/pkgs/gridBase_0.4-7.tar.gz \
 /opt/pkgs/GSEABase_1.56.0.tar.gz \
 /opt/pkgs/hardhat_1.3.0.tar.gz \
 /opt/pkgs/RhpcBLASctl_0.23-42.tar.gz \
 /opt/pkgs/hdf5r_1.3.5.tar.gz \
 /opt/pkgs/hexbin_1.28.3.tar.gz \
 /opt/pkgs/interactiveDisplayBase_1.32.0.tar.gz \
 /opt/pkgs/survival_3.2-13.tar.gz \
 /opt/pkgs/MASS_7.3.55-tar.gz \
 /opt/pkgs/SQUAREM_2021.1.tar.gz \
 /opt/pkgs/lava_1.7.2.1.tar.gz \
 /opt/pkgs/lobstr_1.1.2.tar.gz \
 /opt/pkgs/lsei_1.3-0.tar.gz \
 /opt/pkgs/markdown_1.13.tar.gz \
 /opt/pkgs/MAST_1.20.0.tar.gz 

# /opt/pkgs/zip_2.3.3.tar.gz \

RUN set -eux; \
    R CMD INSTALL -l /usr/local/lib/R/site-library \
 /opt/pkgs/lpSolve_5.6.16.tar.gz \
 /opt/pkgs/mclust_6.0.0.tar.gz \
 /opt/pkgs/metapod_1.2.0.tar.gz \
 /opt/pkgs/npsurv_0.5-0.tar.gz \
 /opt/pkgs/pheatmap_1.0.12.tar.gz \
 /opt/pkgs/polynom_1.4-1.tar.gz \
 /opt/pkgs/proj4_1.0-12.tar.gz \
 /opt/pkgs/pryr_0.1.5.tar.gz \
 /opt/pkgs/ps_1.8.1.tar.gz \
 /opt/pkgs/RcppParallel_5.1.6.tar.gz \
 /opt/pkgs/Rhdf5lib_1.16.0.tar.gz \
 /opt/pkgs/rstatix_0.7.0.tar.gz \
 /opt/pkgs/scuttle_1.4.0.tar.gz

RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/scater_1.22.0.tar.gz
RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/statmod_1.5.0.tar.gz
RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/scran_1.22.1.tar.gz
RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/SingleR_1.8.1.tar.gz
RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/svglite_2.1.0.tar.gz
RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/svMisc_1.2.3.tar.gz
RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/sys_3.4.2.tar.gz
RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/timeDate_4022.108.tar.gz
RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/TrajectoryUtils_1.2.0.tar.gz
RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/TSCAN_1.32.0.tar.gz
RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/xts_0.12.2.tar.gz
RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/TTR_0.24.3.tar.gz
RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/xml2_1.3.3.tar.gz

RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/maps_3.4.1.tar.gz
RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/ggpubr_0.4.0.tar.gz
RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/harmony_0.1.1.tar.gz
RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/prodlim_2023.03.31.tar.gz
RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/ipred_0.9-14.tar.gz
RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/recipes_1.0.6.tar.gz
RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/caret_6.0-94.tar.gz
RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/quantmod_0.4.20.tar.gz
RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/rhdf5filters_1.6.0.tar.gz
RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/rhdf5_2.38.1.tar.gz
RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/HDF5Array_1.22.1.tar.gz
RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/GSVA_1.42.0.tar.gz
RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/AnnotationHub_3.2.2.tar.gz
RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/ExperimentHub_2.2.1.tar.gz
RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/celldex_1.4.0.tar.gz

RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/xgboost_1.7.8.1.tar.gz
RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/scDblFinder_1.8.0.tar.gz
# new
RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/profvis_0.3.7.tar.gz
RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/sessioninfo_1.2.3.tar.gz

RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/downlit_0.4.2.tar.gz
RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/whisker_0.4.1.tar.gz
RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/pkgdown_2.0.7.tar.gz

RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/pkgbuild_1.4.8.tar.gz


RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/brew_1.0-8.tar.gz
RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/roxygen2_7.2.3.tar.gz
RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/urlchecker_1.0.1.tar.gz
RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/rversions_2.1.2.tar.gz
RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/xopen_1.0.0.tar.gz
RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/rcmdcheck_1.4.0.tar.gz


RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/credentials_1.3.2.tar.gz
RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/gert_1.9.2.tar.gz
RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/gitcreds_0.1.2.tar.gz
RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/httr2_0.2.3.tar.gz
RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/ini_0.3.1.tar.gz

RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/gh_1.4.0.tar.gz

RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/usethis_3.1.0.tar.gz
RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/devtools_2.4.5.tar.gz

#deprecated
RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/ggalt_0.4.0.tar.gz
RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/EnhancedVolcano_1.12.0.tar.gz
RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/nidapFunctions_0.7.8.tar.gz

# RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/SparseArray_1.2.0.tar.gz

# _______ end 

# hold off for now 
#RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/SCWorkflow_NA.tar.gz

# hard  - difficult dependencies
#RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/litedown_0.1.tar.gz
#RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/litedown_0.7.tar.gz
#RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/SparseArray_1.8.0.tar.gz
# genomation ?  what for? it's just data for genomation
#RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/genomationData_1.24.0.tar.gz

# not quite version ... use older version 
#RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/xgboost_1.7.8.1.tar.gz

# RUN R -e "install.packages('https://github.com/NIDAP-Community/nidapFunctions/raw/main/nidapFunctions_0.7.8.tar.gz', repos=NULL)"

# problems:
# RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/vector_0.0.2.tar.gz   NO 
# RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/arrow_4.0.1.tar.gz   NO 
# RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/FoundryObjects_0.86.0.tar.gz   NO 
# RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/FoundrySparkR_0.2.0.tar.gz   NO 
# RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/SparkR_3.4.1.tar.gz   NO 

#dupes:
# /opt/pkgs/viridis_0.6.5.tar.gz \
# RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/plotly_4.10.4.tar.gz
# RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/dplyr_1.1.4.tar.gz
#RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/GenomeInfoDbData_1.2.7.tar.gz

#RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/SCWorkflow_NA.tar.gz

RUN R CMD INSTALL -l /usr/local/lib/R/site-library /opt/pkgs/remotes_2.5.0.tar.gz
#        RUN R -e "remotes::install_github('NIDAP-Community/SCWorkflow/tree/GalaxyCLI')"
# worked RUN R -e "remotes::install_github('https://github.com/NIDAP-Community/SCWorkflow/tree/GalaxyCLI')"

RUN R -e "remotes::install_github('NIDAP-Community/SCWorkflow', ref = 'GalaxyCLI')"

COPY Dockerfile /

RUN rm -rf /opt/pkgs


