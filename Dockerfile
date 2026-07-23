FROM bioconductor/bioconductor_docker:latest
MAINTAINER Disease Transcriptomics Lab <diseasetranscriptomicslab@gmail.com>

RUN apt-get update && apt-get -y upgrade && apt-get -y autoremove

RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages('remotes')"

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

# # To start a plain R/RStudio session instead (e.g. for debugging):
# docker run -ti [docker image] R
# docker run -e PASSWORD=bioc -p 8787:8787 [docker image]   # RStudio on :8787
