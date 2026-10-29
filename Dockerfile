# processor-flowsom: FlowSOM clustering of transformed FCS files (R, Bioconductor).
# R 4.6.0 (rocker; CRAN pinned to the Posit Public Package Manager snapshot of
# 2026-06-23 by the base image) + Bioconductor 3.23 binaries.
FROM rocker/r-ver:4.6.0@sha256:6f05a1a8b8c52328f181593923909b01cbfd14c9caea93bf75ddc65e806d8eac

ENV BIOC_VERSION=3.23 \
    DEBIAN_FRONTEND=noninteractive

# Shared libraries the binary packages link against (igraph: GLPK/GMP/libxml2;
# XML: libxml2; flowCore: libcurl-gnutls), plus a font for the PNG figures.
RUN apt-get update \
 && apt-get install -y --no-install-recommends libxml2 libglpk40 libgmp10 libcurl3t64-gnutls fonts-dejavu-core \
 && rm -rf /var/lib/apt/lists/*

# Bioconductor 3.23 packages come from Bioconductor's prebuilt binaries for its
# rocker-based containers; CRAN packages from the pinned P3M snapshot. The build
# fails if a version differs from what the README and provenance document.
RUN Rscript -e ' \
  bioc <- paste0("https://bioconductor.org/packages/", Sys.getenv("BIOC_VERSION")); \
  options(repos = c(BiocBinaries = paste0(bioc, "/container-binaries/bioconductor_docker"), \
                    BioCsoft = paste0(bioc, "/bioc"), \
                    BioCexp = paste0(bioc, "/data/experiment"), \
                    BioCann = paste0(bioc, "/data/annotation"), \
                    CRAN = Sys.getenv("CRAN")), Ncpus = 4); \
  install.packages(c("FlowSOM", "jsonlite")); \
  for (p in c("FlowSOM", "flowCore", "jsonlite", "ggplot2", "igraph", "ConsensusClusterPlus")) \
    suppressPackageStartupMessages(library(p, character.only = TRUE)); \
  stopifnot(packageVersion("FlowSOM") == "2.20.0", packageVersion("flowCore") == "2.24.0"); \
  ip <- installed.packages(fields = "License")[, c("Package", "Version", "License")]; \
  write.csv(ip, "/usr/local/lib/R/installed-packages.csv", row.names = FALSE)' \
 && rm -rf /tmp/Rtmp* /tmp/downloaded_packages \
 && cd /usr/local/lib/R/site-library && rm -rf */doc */html */help */tests */examples

WORKDIR /app
COPY R/ /app/R/
COPY VERSION /app/VERSION

# Runs as uid 1000 on the platform; nothing here needs root.
USER 1000
ENTRYPOINT ["Rscript", "/app/R/main.R"]
