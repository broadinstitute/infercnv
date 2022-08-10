# Docker file for inferCNV
FROM bioconductor/bioconductor_docker:devel

LABEL org.label-schema.license="BSD-3-Clause" \
      org.label-schema.vendor="Broad Institute" \
      maintainer="Christophe Georgescu <cgeorges@broadinstitute.org>"

RUN apt-get update && apt-get -y install curl libssl-dev libcurl4-openssl-dev \
                                        libxml2-dev git python3 jags \
                                        r-cran-rjags time && \
                      apt-get clean && rm -rf /var/tmp/* \
                                          /tmp/* /var/lib/apt/lists/*

# Install R and Bioconductor packages
RUN echo "options(repos = c(CRAN = 'https://cran.rstudio.com'))" >.Rprofile
RUN R -e "BiocManager::install('infercnv')"
#RUN R -e "install.packages(c('cluster', 'Seurat', 'parallelDist', 'optparse'), repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages(c('cluster', 'Seurat', 'optparse', 'igraph', 'reshape2'), repos = 'http://cran.us.r-project.org')"
#RUN R -e "install.packages('phyclust', repos = 'http://cran.us.r-project.org')"

# RUN pip3 install numpy igraph pandas leidenalg
# ENV RETICULATE_PYTHON=/usr/bin/python3

# Checkout and install infercnv
# update to version bump commit
RUN git clone https://github.com/broadinstitute/infercnv && cd infercnv && \
      git checkout master && git checkout 1b46b48303bac4a882bcb758e78fcf7f832fdefb && \
      R CMD INSTALL . 

ENV PATH=${PATH}:/infercnv/scripts

CMD inferCNV.R --help

