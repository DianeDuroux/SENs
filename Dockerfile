FROM rocker/r-ver:4.3.1

# Install system dependencies needed for R packages
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    libgsl-dev \
    liblapack-dev \
    libblas-dev \
    gfortran \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Copy repo into container
COPY . /app

# Install CRAN packages
RUN R -e "install.packages(c( \
    'GeneCycle', \
    'huge', \
    'combinat', \
    'NetworkDistance', \
    'graphkernels', \
    'progress', \
    'stringr', \
    'tidyr', \
    'KRLS', \
    'anocva', \
    'data.table', \
    'readr', \
    'ggplot2', \
    'dplyr', \
    'reshape2', \
    'RColorBrewer', \
    'ggdendro', \
    'Matrix', \
    'SparseM', \
    'pracma', \
    'poLCA' \
    ), repos = 'https://cloud.r-project.org')"

# Install Bioconductor packages (org.Hs.eg.db + AnnotationDbi)
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org'); \
          BiocManager::install(c('org.Hs.eg.db', 'AnnotationDbi'), ask = FALSE)"

# Make sure output folders exist
RUN mkdir -p results/significance/Gene results/netANOVA results/plots results/aggregated_network

# Default run command (can be overridden)
CMD ["Rscript", "run_pipeline.R"]
