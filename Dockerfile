FROM rocker/r-ver:latest
MAINTAINER Disease Transcriptomics Lab <diseasetranscriptomicslab@gmail.com>

# rocker/r-ver is the same lineage as bioconductor/bioconductor_docker
# (which this Dockerfile used before) but WITHOUT RStudio Server bolted on
# top - this app is a headless Shiny app and never uses it, so skipping it
# saves close to 1.7GB. The trade-off: bioconductor_docker also came with
# BiocManager and various system libraries pre-installed for building
# Bioconductor/CRAN packages from source; both are added explicitly below
# since rocker/r-ver doesn't include them.
RUN apt-get update && apt-get -y upgrade && apt-get install -y --no-install-recommends \
      build-essential \
      gfortran \
      pkg-config \
      libcurl4-openssl-dev \
      libssl-dev \
      libxml2-dev \
      libgit2-dev \
      libuv1-dev \
      zlib1g-dev \
      libpng-dev \
      liblapack-dev \
      libblas-dev \
      libglpk-dev \
    && apt-get -y autoremove

RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages(c('remotes', 'BiocManager'))"

# Copy package source code
WORKDIR /markeR
ADD . .

# Install markeR and all of its Imports (Bioconductor + CRAN) from source
RUN Rscript -e "remotes::install_local()"

# markeR's gene-ID-conversion feature (Preprocessing tab) uses biomaRt
# (works for all species via the Ensembl API) with AnnotationDbi + a local
# org.*.eg.db package as an offline fallback. Both are only in DESCRIPTION's
# Suggests (not Imports), so remotes::install_local() above does not
# reliably pull them in - install them explicitly so the feature works out
# of the box in the deployed app instead of just warning users to install
# something they have no way to install themselves.
RUN Rscript -e "BiocManager::install(c('biomaRt', 'AnnotationDbi', 'org.Hs.eg.db', 'org.Mm.eg.db'), update = FALSE, ask = FALSE)"

EXPOSE 3838

# Launch the markeR Shiny app when the container starts, listening on all
# interfaces so it's reachable from outside the container.
CMD ["R", "-e", "markeR::markeRapp(host = '0.0.0.0', port = 3838, launch.browser = FALSE)"]

# # To start a plain R session instead (e.g. for debugging):
# docker run -ti [docker image] R
