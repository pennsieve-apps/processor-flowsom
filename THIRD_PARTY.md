# Third-party software

Our code (`R/`, `tests/`) is Apache-2.0 (see `LICENSE`). The image runs it with R and
R packages that are GPL-licensed and loaded into the same process (FlowSOM is
GPL (>= 2); `cytolib`, which flowCore links, is AGPL-3.0-only), so **the image as a whole
is distributed under the GPL/AGPL terms of those components**; each package keeps its own
licence, and their sources are on CRAN / Bioconductor at the versions below. Nothing here
is modified.

| Component | Version | What | Licence | Source |
|---|---|---|---|---|
| FlowSOM | 2.20.0 | SOM clustering, metaclustering, mapping (`NewData`), star tree | GPL (>= 2) | https://bioconductor.org/packages/FlowSOM (Bioconductor 3.23) |
| flowCore | 2.24.0 | FCS reading and writing | Artistic-2.0 | https://bioconductor.org/packages/flowCore |
| cytolib | 2.24.0 | C++ cytometry library under flowCore | AGPL-3.0-only | https://bioconductor.org/packages/cytolib |
| ConsensusClusterPlus | 1.76.0 | Metaclustering (called by FlowSOM) | GPL-2 | https://bioconductor.org/packages/ConsensusClusterPlus |
| igraph | 2.3.3 | Minimal spanning tree (called by FlowSOM) | GPL (>= 2) | CRAN |
| ggplot2 | 4.0.3 | Heatmap, star tree | MIT | CRAN |
| jsonlite | 2.0.0 | JSON in/out | MIT | CRAN |
| R | 4.6.0 | Interpreter (`rocker/r-ver:4.6.0`, Ubuntu 24.04) | GPL-2 \| GPL-3 | https://www.r-project.org, https://rocker-project.org |
| Ubuntu packages | noble | libxml2, GLPK, GMP, libcurl-gnutls, DejaVu fonts | their own free licences (MIT, GPL-3, LGPL-3, curl, Bitstream Vera) | https://packages.ubuntu.com |

CRAN packages come from the Posit Public Package Manager snapshot of 2026-06-23 (pinned by
the base image); Bioconductor 3.23 packages from Bioconductor's prebuilt container binaries.
The image lists every installed package with its version and licence in
`/usr/local/lib/R/installed-packages.csv`. All of them (R's base and recommended
packages aside):

abind 1.4-8 (MIT + file LICENSE); ALL 1.54.0 (Artistic-2.0); backports 1.5.1 (GPL-2); base64enc 0.1-6 (GPL-2); BH 1.90.0-1 (BSL-1.0); Biobase 2.72.0 (Artistic-2.0); BiocGenerics 0.58.1 (Artistic-2.0); biocmake 1.4.0 (MIT + file LICENSE); boot 1.3-32 (Unlimited); broom 1.0.13 (MIT + file LICENSE); car 3.1-5 (GPL (>= 2)); carData 3.0-6 (GPL (>= 2)); class 7.3-23 (GPL-2); cli 3.6.6 (MIT + file LICENSE); cluster 2.1.8.2 (GPL (>= 2)); codetools 0.2-20 (GPL); colorRamps 2.3.4 (GPL); colorspace 2.1-3 (BSD_3_clause + file LICENSE); ConsensusClusterPlus 1.76.0 (GPL version 2); corrplot 0.95 (MIT + file LICENSE); cowplot 1.2.0 (GPL-2); cpp11 0.5.5 (MIT + file LICENSE); cytolib 2.24.0 (AGPL-3.0-only); Deriv 4.2.0 (GPL (>= 3)); dir.expiry 1.20.0 (GPL-3); doBy 4.7.2 (GPL (>= 2)); docopt 0.7.2 (MIT + file LICENSE); dplyr 1.2.1 (MIT + file LICENSE); farver 2.1.2 (MIT + file LICENSE); filelock 1.0.3 (MIT + file LICENSE); flowCore 2.24.0 (Artistic-2.0); FlowSOM 2.20.0 (GPL (>= 2)); forecast 9.0.2 (GPL-3); foreign 0.8-91 (GPL (>= 2)); Formula 1.2-5 (GPL-2); fracdiff 1.5-4 (GPL (>= 2)); generics 0.1.4 (MIT + file LICENSE); ggforce 0.5.0 (MIT + file LICENSE); ggnewscale 0.5.2 (GPL-3); ggplot2 4.0.3 (MIT + file LICENSE); ggpubr 1.0.0 (GPL (>= 2)); ggrepel 0.9.8 (GPL-3); ggsci 5.1.0 (GPL (>= 3)); ggsignif 0.6.4 (GPL-3); glue 1.8.1 (MIT + file LICENSE); gridExtra 2.3.1 (GPL (>= 2)); gtable 0.3.6 (MIT + file LICENSE); igraph 2.3.3 (GPL (>= 2)); isoband 0.3.0 (MIT + file LICENSE); jsonlite 2.0.0 (MIT + file LICENSE); KernSmooth 2.23-26 (Unlimited); labeling 0.4.3 (MIT + file LICENSE); lattice 0.22-9 (GPL (>= 2)); lifecycle 1.0.5 (MIT + file LICENSE); littler 0.3.23 (GPL (>= 2)); lme4 2.0-6 (GPL (>= 2)); lmtest 0.9-40 (GPL-2); magrittr 2.0.5 (MIT + file LICENSE); MASS 7.3-65 (GPL-2); Matrix 1.7-5 (GPL (>= 2)); MatrixModels 0.5-4 (GPL (>= 2)); matrixStats 1.5.0 (Artistic-2.0); mgcv 1.9-4 (GPL (>= 2)); minqa 1.2.8 (GPL-2); modelr 0.1.11 (GPL-3); nlme 3.1-169 (GPL (>= 2)); nloptr 2.2.1 (LGPL (>= 3)); nnet 7.3-20 (GPL-2); numDeriv 2016.8-1.1 (GPL-2); pbkrtest 0.5.5 (GPL (>= 2)); pillar 1.11.1 (MIT + file LICENSE); pkgconfig 2.0.3 (MIT + file LICENSE); polyclip 1.10-7 (BSL); polynom 1.4-1 (GPL-2); purrr 1.2.2 (MIT + file LICENSE); quantreg 6.1 (GPL (>= 2)); R6 2.6.1 (MIT + file LICENSE); rbibutils 2.4.1 (GPL-2); RColorBrewer 1.1-3 (Apache License 2.0); Rcpp 1.1.1-1.1 (GPL (>= 2)); RcppArmadillo 15.4.0-1 (GPL (>= 2)); RcppEigen 0.3.4.0.2 (GPL (>= 2)); Rdpack 2.6.6 (GPL (>= 2)); reformulas 0.4.4 (GPL-3); Rhdf5lib 2.0.0 (Artistic-2.0); rlang 1.2.0 (MIT + file LICENSE); rpart 4.1.27 (GPL-2); RProtoBufLib 2.24.0 (BSD_3_clause); rstatix 1.0.0 (GPL-2); Rtsne 0.17 (file LICENSE); S4Vectors 0.50.3 (Artistic-2.0); S7 0.2.2 (MIT + file LICENSE); scales 1.4.0 (MIT + file LICENSE); SparseM 1.84-2 (GPL (>= 2)); spatial 7.3-18 (GPL-2); stringi 1.8.7 (file LICENSE); stringr 1.6.0 (MIT + file LICENSE); survival 3.8-6 (LGPL (>= 2)); systemfonts 1.3.2 (MIT + file LICENSE); tibble 3.3.1 (MIT + file LICENSE); tidyr 1.3.2 (MIT + file LICENSE); tidyselect 1.2.1 (MIT + file LICENSE); timeDate 4052.112 (GPL (>= 2)); tweenr 2.0.3 (MIT + file LICENSE); urca 1.3-4 (GPL (>= 2)); utf8 1.2.6 (Apache License (== 2.0)); vctrs 0.7.3 (MIT + file LICENSE); viridisLite 0.4.3 (MIT + file LICENSE); withr 3.0.3 (MIT + file LICENSE); XML 3.99-0.23 (BSD_3_clause + file LICENSE); zoo 1.8-15 (GPL-2).

`tests/make_fixture.py` (not in the image) uses NumPy (BSD-3) and FlowIO (BSD-3).

## Citations

See README.md, Citations.
