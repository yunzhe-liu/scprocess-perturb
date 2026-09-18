#!/usr/bin/env bash
set -euo pipefail
Rscript -e "remotes::install_github('immunogenomics/presto@3c97180d90d330524b5b8353a342ab2f4b990d21', upgrade='never', dependencies=FALSE, repos='https://cloud.r-project.org')"
Rscript -e "stopifnot(packageVersion('presto') == '1.1.0')"
