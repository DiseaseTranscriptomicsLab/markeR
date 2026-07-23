FROM bioconductor/bioconductor_docker:latest
MAINTAINER Disease Transcriptomics Lab <diseasetranscriptomicslab@gmail.com>

RUN apt-get update && apt-get -y upgrade && apt-get -y autoremove

RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages('remotes')"

WORKDIR /markeR

# ---- Dependency layer (cached) ----
# Copy ONLY the dependency manifest first and install everything it lists.
# Docker only re-runs this (slow, ~9 min) step if DESCRIPTION itself changes.
# Editing R code or the bundled example data below no longer invalidates it.
# dependencies = NA (not TRUE): installs Depends/Imports/LinkingTo only, not
# the full Suggests tree (devtools, roxygen2, covr, BiocStyle, renv, sva,
# tximport, etc.) - those are for package development/checking, not needed
# to run the deployed Shiny app, and installing them was what pushed the
# build from ~8 min to 16+ min.
COPY DESCRIPTION /markeR/DESCRIPTION
RUN Rscript -e "remotes::install_deps(dependencies = NA)"

# markeR's gene-ID-conversion feature (Preprocessing tab) uses biomaRt
# (works for all species via the Ensembl API) with AnnotationDbi + a local
# org.*.eg.db package as an offline fallback. Both are only in DESCRIPTION's
# Suggests (not Imports), so install_deps() above does not reliably pull
# them in - install them explicitly so the feature works out of the box in
# the deployed app instead of just warning users to install something they
# have no way to install themselves. Also part of the cached dependency
# layer since it doesn't depend on the source code either.
RUN Rscript -e "BiocManager::install(c('biomaRt', 'AnnotationDbi', 'org.Hs.eg.db', 'org.Mm.eg.db'), update = FALSE, ask = FALSE)"

# ---- Source layer (rebuilt whenever code/data changes) ----
# Now copy the actual package source (code + bundled example data under
# inst/appdata). All dependencies are already installed above, so this only
# has to install the local package itself - fast, regardless of source size.
ADD . .
RUN Rscript -e "remotes::install_local()"

EXPOSE 3838

# Launch the markeR Shiny app when the container starts, listening on all
# interfaces so it's reachable from outside the container.
CMD ["R", "-e", "markeR::markeRapp(host = '0.0.0.0', port = 3838, launch.browser = FALSE)"]

# # To start a plain R/RStudio session instead (e.g. for debugging):
# docker run -ti [docker image] R
# docker run -e PASSWORD=bioc -p 8787:8787 [docker image]   # RStudio on :8787
