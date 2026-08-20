FROM bioconductor/bioconductor_docker:latest
MAINTAINER Disease Transcriptomics Lab <diseasetranscriptomicslab@gmail.com>

RUN apt-get update && apt-get -y upgrade && apt-get -y autoremove && apt-get clean

RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages('remotes')"

WORKDIR /markeR

# Copy ONLY the dependency manifest first, before the rest of the source
# code. Docker reuses a layer as long as its inputs haven't changed, so as
# long as DESCRIPTION doesn't change, this layer (and the two install steps
# below) are reused even when you edit R/Shiny source files - meaning a
# two-line script change no longer triggers a full reinstall of every
# Bioconductor/CRAN dependency from scratch.
COPY DESCRIPTION .

RUN Rscript -e "remotes::install_deps(dependencies = c('Depends', 'Imports', 'LinkingTo'))"

# A handful of packages are only in DESCRIPTION's Suggests (not Imports),
# so remotes::install_deps() above does not reliably pull them in - but
# they're genuinely used at runtime (each guarded behind
# requireNamespace(..., quietly = TRUE) so the app degrades gracefully if
# they're missing, rather than failing outright). Install them explicitly
# so these features work out of the box in the deployed app instead of
# just telling users to install something they have no way to install
# themselves:
#   - biomaRt (+ AnnotationDbi, org.Hs.eg.db, org.Mm.eg.db as an offline
#     fallback) powers the gene-ID-conversion feature (Preprocessing tab).
#   - stringdist powers fuzzy sample-ID matching in the GEO import module.
#   - tximport is required to import Salmon .sf files in the GEO import
#     module; without it that import path hard-fails.
RUN Rscript -e "BiocManager::install(c('biomaRt', 'AnnotationDbi', 'org.Hs.eg.db', 'org.Mm.eg.db', 'stringdist', 'tximport'), update = FALSE, ask = FALSE)"

# NOW copy the rest of the source code and install markeR itself. This step
# no longer needs to resolve or download any dependencies (they're already
# installed above), so it's fast and cheap - ordinary code edits only redo
# this quick final step, not the whole dependency tree.
ADD . .
RUN Rscript -e "remotes::install_local(dependencies = FALSE)"

EXPOSE 3838

# Launch the markeR Shiny app when the container starts, listening on all
# interfaces so it's reachable from outside the container.
CMD ["R", "-e", "markeR::markeRapp(host = '0.0.0.0', port = 3838, launch.browser = FALSE)"]

# # To start a plain R/RStudio session instead (e.g. for debugging):
# docker run -ti [docker image] R
# docker run -e PASSWORD=bioc -p 8787:8787 [docker image]   # RStudio on :8787
