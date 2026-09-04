#!/usr/bin/env bash

# Bioconda stops at xTail 1.1.5, which fails with the R-4-era
# SummarizedExperiment API. Install the official upstream compatibility release
# into this already activated, Snakemake-managed environment.
set -euo pipefail

readonly source_url="https://github.com/xryanglab/xtail/releases/download/v1.2.0/xtail_1.2.0-source.tar.gz"
readonly source_sha256="5975e7b9ea692be69ebaee85a64f710ca742acf6f2464d961c764fe5f1cccc95"
: "${CONDA_PREFIX:?xTail post-deploy requires an activated Conda environment}"
readonly target_library="$CONDA_PREFIX/lib/R/library"

deploy_tmp=$(mktemp -d "${TMPDIR:-/tmp}/hribo-xtail.XXXXXX")
cleanup() {
    rm -rf -- "$deploy_tmp"
}
trap cleanup EXIT

archive="$deploy_tmp/xtail_1.2.0-source.tar.gz"
curl --fail --location --retry 3 --proto '=https' --tlsv1.2 \
    --output "$archive" "$source_url"
printf '%s  %s\n' "$source_sha256" "$archive" | sha256sum --check -

R --vanilla CMD INSTALL --clean --no-multiarch --library="$target_library" "$archive"
Rscript --vanilla -e '
library_path <- normalizePath(commandArgs(trailingOnly = TRUE)[1], mustWork = TRUE)
package_path <- normalizePath(find.package("xtail", lib.loc = library_path), mustWork = TRUE)
stopifnot(
  as.character(packageVersion("xtail", lib.loc = library_path)) == "1.2.0",
  identical(dirname(package_path), library_path)
)
' "$target_library"
