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

# Install packages with specific versions from CRAN and Bioconductor
RUN R --vanilla --slave -e " \
  options(repos = c(CRAN = 'https://cran.r-project.org')); \
  \
  if (!require('BiocManager', quietly = TRUE)) { \
    install.packages('BiocManager', quiet = TRUE) \
  }; \
  if (!require('remotes', quietly = TRUE)) { \
    install.packages('remotes', quiet = TRUE) \
  }; \
  \
  # Install specific CRAN package versions \
  cran_versions <- list( \
    bitops = '1.0-8', \
    RCurl = '1.98-1.14', \
    GlobalOptions = '0.1.2', \
    bit = '4.0.4', \
    assertthat = '0.2.1', \
    bit64 = '4.0.5', \
    shape = '1.4.6', \
    colorspace = '2.0-3', \
    circlize = '0.4.15', \
    rlang = '1.1.4', \
    cli = '3.6.3', \
    glue = '1.7.0', \
    lifecycle = '1.0.4', \
    vctrs = '0.6.5', \
    blob = '1.2.3', \
    clue = '0.3-61', \
    cluster = '2.1.2', \
    codetools = '0.2-18', \
    crayon = '1.5.3', \
    data.table = '1.15.4', \
    DBI = '1.2.3', \
    dendsort = '0.3.4', \
    digest = '0.6.37', \
    ellipsis = '0.3.2', \
    evaluate = '0.24.0', \
    fansi = '1.0.6', \
    farver = '2.1.1', \
    fastmap = '1.2.0', \
    cachem = '1.1.0', \
    fastmatch = '1.1-3', \
    iterators = '1.0.14', \
    foreach = '1.5.2', \
    formatR = '1.14', \
    generics = '0.1.3', \
    rjson = '0.2.21', \
    GetoptLong = '1.0.5', \
    gtable = '0.3.5', \
    labeling = '0.4.2', \
    munsell = '0.5.0', \
    R6 = '2.5.1', \
    gridExtra = '2.3', \
    RColorBrewer = '1.1-3', \
    utf8 = '1.2.4', \
    pillar = '1.9.0', \
    pkgconfig = '2.0.3', \
    viridisLite = '0.4.1', \
    scales = '1.2.1', \
    withr = '3.0.1', \
    gtools = '3.9.5', \
    gridGraphics = '0.5-1', \
    hms = '1.1.2', \
    irlba = '2.3.5.1', \
    jsonlite = '1.8.8', \
    png = '0.1-7', \
    KernSmooth = '2.23-20', \
    lambda.r = '1.2.4', \
    lattice = '0.20-45', \
    lazyeval = '0.2.2', \
    locfit = '1.5-9.9', \
    Matrix = '1.5-1', \
    matrixStats = '0.62.0', \
    memoise = '2.0.1', \
    mgcv = '1.8-39', \
    nlme = '3.1-155', \
    rsvd = '1.0.5', \
    Rcpp = '1.0.13', \
    uuid = '1.1-0', \
    xfun = '0.47', \
    xtable = '1.8-4', \
    yaml = '2.3.10', \
    beeswarm = '0.4.0', \
    filelock = '1.0.3', \
    bslib = '0.4.2', \
    Cairo = '1.6-0', \
    desc = '1.4.2', \
    pkgload = '1.3.0', \
    brio = '1.1.3', \
    zip = '2.3.3', \
    openxlsx = '4.2.5.2', \
    praise = '1.0.0', \
    diffobj = '0.3.5', \
    waldo = '0.5.1', \
    testthat = '3.1.6', \
    nloptr = '2.0.3', \
    minqa = '1.2.5', \
    lme4 = '1.1-33', \
    MatrixModels = '0.5-1', \
    SparseM = '1.81', \
    quantreg = '5.95', \
    numDeriv = '2016.8-1.1', \
    pbkrtest = '0.5.2', \
    maptools = '1.1-7', \
    carData = '3.0-5', \
    rio = '0.5.29', \
    car = '3.0-0', \
    pROC = '1.18.2', \
    proxy = '0.4-27', \
    e1071 = '1.7-13', \
    ModelMetrics = '1.2.2.2', \
    clock = '0.7.0', \
    shinyjs = '2.1.0', \
    colourpicker = '1.2.0', \
    combinat = '0.0-8', \
    corrplot = '0.92', \
    cpp11 = '0.4.7', \
    diagram = '1.6.5', \
    DT = '0.28', \
    Rttf2pt1 = '1.3.12', \
    extrafontdb = '1.0', \
    extrafont = '0.18', \
    vipor = '0.4.7', \
    ggbeeswarm = '0.7.2', \
    ggrastr = '1.0.2', \
    fastICA = '1.2-3', \
    gdata = '2.18.0.1', \
    ggExtra = '0.10.1', \
    ggsci = '3.0.0', \
    ggsignif = '0.6.3', \
    gower = '1.0.1', \
    gridBase = '0.4-7', \
    hardhat = '1.3.0', \
    RhpcBLASctl = '0.23-42', \
    hdf5r = '1.3.5', \
    hexbin = '1.28.3', \
    survival = '3.2-13', \
    MASS = '7.3-55', \
    SQUAREM = '2021.1', \
    lava = '1.7.2.1', \
    lobstr = '1.1.2', \
    lsei = '1.3-0', \
    markdown = '1.13', \
    lpSolve = '5.6.16', \
    mclust = '6.0.0', \
    npsurv = '0.5-0', \
    pheatmap = '1.0.12', \
    polynom = '1.4-1', \
    proj4 = '1.0-12', \
    pryr = '0.1.5', \
    ps = '1.8.1', \
    RcppParallel = '5.1.6', \
    rstatix = '0.7.0', \
    timeDate = '4022.108', \
    maps = '3.4.1', \
    ggpubr = '0.4.0', \
    prodlim = '2023.03.31', \
    ipred = '0.9-14', \
    recipes = '1.0.6', \
    caret = '6.0-94', \
    quantmod = '0.4.20', \
    profvis = '0.3.7', \
    sessioninfo = '1.2.3', \
    downlit = '0.4.2', \
    whisker = '0.4.1', \
    pkgdown = '2.0.7', \
    pkgbuild = '1.4.8', \
    brew = '1.0-8', \
    roxygen2 = '7.2.3', \
    urlchecker = '1.0.1', \
    rversions = '2.1.2', \
    xopen = '1.0.0', \
    rcmdcheck = '1.4.0', \
    credentials = '1.3.2', \
    gert = '1.9.2', \
    gitcreds = '0.1.2', \
    httr2 = '0.2.3', \
    ini = '0.3.1', \
    gh = '1.4.0', \
    usethis = '3.1.0', \
    devtools = '2.4.5', \
    ggalt = '0.4.0', \
    EnhancedVolcano = '1.12.0' \
  ); \
  \
  for (pkg in names(cran_versions)) { \
    tryCatch({ \
      remotes::install_version(pkg, version = cran_versions[[pkg]], repos = 'https://cran.r-project.org', quiet = TRUE) \
    }, error = function(e) { \
      message('Note: install_version failed for ', pkg, ', trying install.packages') \
      install.packages(pkg, quiet = TRUE) \
    }) \
  }; \
  \
  # Install Bioconductor packages with specific versions \
  bioc_versions <- list( \
    BiocGenerics = '0.40.0', \
    Biobase = '2.54.0', \
    S4Vectors = '0.32.4', \
    IRanges = '2.28.0', \
    XVector = '0.34.0', \
    GenomeInfoDbData = '1.2.7', \
    GenomeInfoDb = '1.30.1', \
    beachmat = '2.10.0', \
    edgeR = '3.36.0', \
    GenomicRanges = '1.46.1', \
    Biostrings = '2.62.0', \
    DelayedArray = '0.20.0', \
    sparseMatrixStats = '1.6.0', \
    DelayedMatrixStats = '1.16.0', \
    ScaledMatrix = '1.2.0', \
    zlibbioc = '1.40.0', \
    KEGGREST = '1.34.0', \
    AnnotationDbi = '1.56.2', \
    BiocParallel = '1.28.3', \
    BiocSingular = '1.10.0', \
    ComplexHeatmap = '2.10.0', \
    fgsea = '1.20.0', \
    BiocNeighbors = '1.12.0', \
    BiocFileCache = '2.2.1', \
    bluster = '1.4.0', \
    scuttle = '1.4.0', \
    scater = '1.22.0', \
    scran = '1.22.1', \
    SingleR = '1.8.1', \
    MAST = '1.20.0', \
    scDblFinder = '1.8.0', \
    SummarizedExperiment = '1.24.0', \
    SingleCellExperiment = '1.16.0', \
    HDF5Array = '1.22.1', \
    rhdf5 = '2.38.1', \
    rhdf5filters = '1.6.0', \
    GSVA = '1.42.0', \
    ExperimentHub = '2.2.1', \
    AnnotationHub = '3.2.2', \
    celldex = '1.4.0', \
    annotate = '1.72.0', \
    graph = '1.72.0', \
    GSEABase = '1.56.0', \
    interactiveDisplayBase = '1.32.0', \
    TrajectoryUtils = '1.2.0', \
    TSCAN = '1.32.0', \
    conquer = '1.3.3', \
    metapod = '1.2.0', \
    statmod = '1.5.0' \
  ); \
  \
  for (pkg in names(bioc_versions)) { \
    tryCatch({ \
      remotes::install_version(pkg, version = bioc_versions[[pkg]], repos = BiocManager::repositories(), quiet = TRUE) \
    }, error = function(e) { \
      message('Note: install_version failed for ', pkg, ', trying BiocManager') \
      BiocManager::install(pkg, ask = FALSE, quiet = TRUE) \
    }) \
  }; \
  \
  message('All packages installed successfully')"


# Additional packages that need special installation

# Install SCWorkflow from GitHub
RUN R --vanilla --slave -e "remotes::install_github('NIDAP-Community/SCWorkflow', ref = 'GalaxyCLI', quiet = TRUE)"

# Install spatstat family packages with specific versions
RUN R --vanilla --slave -e " \
  options(repos = c(CRAN = 'https://cran.r-project.org')); \
  if (!require('remotes', quietly = TRUE)) install.packages('remotes', quiet = TRUE); \
  spatstat_versions <- list( \
    spatstat.utils = '3.1-0', \
    spatstat.data = '3.0-0', \
    deldir = '2.0-4', \
    polyclip = '1.10-7', \
    spatstat.univar = '2.0-3', \
    spatstat.geom = '3.0-3', \
    spatstat.random = '2.2-0', \
    abind = '1.4-5', \
    tensor = '1.5', \
    goftest = '1.2-3', \
    spatstat.sparse = '3.1-0', \
    spatstat.core = '2.4-2', \
    spatstat.linnet = '2.2-1' \
  ); \
  for (pkg in names(spatstat_versions)) { \
    tryCatch({ \
      remotes::install_version(pkg, version = spatstat_versions[[pkg]], repos = 'https://cran.r-project.org', quiet = TRUE) \
    }, error = function(e) { \
      message('Note: install_version for ', pkg, ' failed, trying install.packages') \
      install.packages(pkg, quiet = TRUE) \
    }) \
  }"

# Install remaining specialized packages
RUN R --vanilla --slave -e " \
  options(repos = c(CRAN = 'https://cran.r-project.org')); \
  if (!require('remotes', quietly = TRUE)) install.packages('remotes', quiet = TRUE); \
  special_versions <- list( \
    RSpectra = '0.16-1', \
    dotCall64 = '1.1-1', \
    spam = '2.11-0', \
    RcppHNSW = '0.4.1', \
    leidenbase = '0.1.30', \
    fastDummies = '1.7.3', \
    sp = '1.5-0', \
    rgeos = '0.5-9', \
    SeuratObject = '4.1.1', \
    Seurat = '4.1.1', \
    XML = '3.99-0.14', \
    anndata = '0.7.5.2', \
    ash = '1.0-15', \
    viridis = '0.6.5', \
    dendextend = '1.16.0', \
    ggrepel = '0.9.5', \
    l2p = '0.0-13', \
    l2psupp = '0.0-13', \
    ica = '1.0-3', \
    Rtsne = '0.16', \
    ggridges = '0.5.3', \
    scattermore = '1.2', \
    listenv = '0.9.1', \
    globals = '0.16.3', \
    parallelly = '1.38.0', \
    future = '1.34.0', \
    future.apply = '1.11.2', \
    RcppEigen = '0.3.3.9.3', \
    RcppAnnoy = '0.0.19', \
    zoo = '1.8-12', \
    lmtest = '0.9-40', \
    fitdistrplus = '1.1-8', \
    caTools = '1.18.2', \
    gplots = '3.1.3', \
    ROCR = '1.0-11', \
    igraph = '2.0.3', \
    pbapply = '1.5-0', \
    commonmark = '1.9.0', \
    httpuv = '1.6.15', \
    sourcetools = '0.1.7-1', \
    shiny = '1.9.1', \
    miniUI = '0.1.1.1', \
    progressr = '0.14.0', \
    sitmo = '2.0.2', \
    dqrng = '0.4.1', \
    FNN = '1.1.3.2', \
    RcppProgress = '0.4.2', \
    uwot = '0.1.14', \
    cowplot = '1.1.1', \
    RcppTOML = '0.2.2', \
    rprojroot = '2.0.4', \
    here = '1.0.1', \
    reticulate = '1.40.0', \
    leiden = '0.4.3', \
    RANN = '2.6.1', \
    RcppArmadillo = '0.12.4.0.0', \
    sctransform = '0.3.4', \
    Rhdf5lib = '1.16.0', \
    xgboost = '1.7.8.1', \
    nidapFunctions = '0.7.8', \
    snow = '0.4-4', \
    BH = '1.81.0-1', \
    futile.options = '1.0.1', \
    futile.logger = '1.4.3', \
    base64enc = '0.1-3', \
    htmltools = '0.5.8.1', \
    rappdirs = '0.3.3', \
    jquerylib = '0.1.4', \
    tinytex = '0.44', \
    fs = '1.6.4', \
    sass = '0.4.6', \
    mime = '0.12', \
    magrittr = '2.0.3', \
    stringi = '1.8.4', \
    stringr = '1.5.1', \
    highr = '0.10', \
    knitr = '1.48', \
    fontawesome = '0.5.1', \
    rmarkdown = '2.28', \
    htmlwidgets = '1.6.4', \
    sys = '3.4.2', \
    askpass = '1.1', \
    openssl = '2.0.5', \
    curl = '6.4.0', \
    httr = '1.4.7', \
    timechange = '0.2.0', \
    lubridate = '1.8.0', \
    isoband = '0.2.7', \
    tibble = '3.2.1', \
    ggplot2 = '3.3.6', \
    patchwork = '1.2.0', \
    later = '1.3.2', \
    promises = '1.3.0', \
    crosstalk = '1.2.0', \
    purrr = '1.0.2', \
    tidyselect = '1.2.1', \
    dplyr = '1.1.4', \
    tidyr = '1.2.1', \
    plotly = '4.10.4', \
    plyr = '1.8.7', \
    clipr = '0.8.0', \
    prettyunits = '1.1.1', \
    progress = '1.2.2', \
    tzdb = '0.4.0', \
    vroom = '1.6.3', \
    readr = '2.1.2', \
    reshape2 = '1.4.4', \
    plogr = '0.2.0', \
    RSQLite = '2.3.9', \
    systemfonts = '1.2.3', \
    textshaping = '0.3.6', \
    ragg = '1.2.5', \
    dbplyr = '2.2.1', \
    rstudioapi = '0.14', \
    dtplyr = '1.3.1', \
    backports = '1.4.1', \
    broom = '1.0.1', \
    processx = '3.8.5', \
    callr = '3.7.3', \
    reprex = '2.0.2', \
    modelr = '0.1.9', \
    conflicted = '1.2.0', \
    rematch2 = '2.1.2', \
    gargle = '1.3.0', \
    rematch = '1.0.1', \
    cellranger = '1.1.0', \
    ids = '1.0.1', \
    googledrive = '2.0.0', \
    googlesheets4 = '1.0.1', \
    readxl = '1.4.1', \
    selectr = '0.4-2', \
    xml2 = '1.3.3', \
    rvest = '1.0.3', \
    forcats = '0.5.2', \
    haven = '2.5.1', \
    tidyverse = '1.3.2', \
    doParallel = '1.0.17', \
    MatrixGenerics = '1.6.0', \
    svglite = '2.1.0', \
    svMisc = '1.2.3', \
    xts = '0.12.2', \
    TTR = '0.24.3' \
  ); \
  for (pkg in names(special_versions)) { \
    tryCatch({ \
      remotes::install_version(pkg, version = special_versions[[pkg]], repos = 'https://cran.r-project.org', quiet = TRUE) \
    }, error = function(e) { \
      message('Note: install_version for ', pkg, ' failed, trying install.packages') \
      install.packages(pkg, quiet = TRUE) \
    }) \
  }"

COPY Dockerfile /


