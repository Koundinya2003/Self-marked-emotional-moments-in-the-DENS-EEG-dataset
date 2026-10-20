#!/usr/bin/env bash
# Builds the R environment used by the pipeline.
# CRAN was not reachable from the original build machine, so R and all packages
# come from conda-forge through micromamba. Quarto 1.5.57 is taken from its
# GitHub release and LaTeX (for the PDF) from the Ubuntu archive.
set -euo pipefail
MM_DIR=${MM_DIR:-/opt/mm}
ENV_DIR=${ENV_DIR:-$MM_DIR/envs/dens}
mkdir -p "$MM_DIR"
if [ ! -x "$MM_DIR/bin/micromamba" ]; then
  curl -sL https://conda.anaconda.org/conda-forge/linux-64/micromamba-2.9.0-0.tar.bz2 | tar -xj -C "$MM_DIR"
fi
export MAMBA_ROOT_PREFIX="$MM_DIR/root"
if [ ! -d "$ENV_DIR" ]; then
  "$MM_DIR/bin/micromamba" create -y -p "$ENV_DIR" -f "$(dirname "$0")/../env/environment.yml"
fi
if [ ! -x /opt/quarto/opt/quarto/bin/quarto ]; then
  curl -sL -o /tmp/quarto.deb https://github.com/quarto-dev/quarto-cli/releases/download/v1.5.57/quarto-1.5.57-linux-amd64.deb
  dpkg -x /tmp/quarto.deb /opt/quarto
fi
if ! command -v xelatex >/dev/null 2>&1; then
  sudo_cmd=""; [ "$(id -u)" -ne 0 ] && sudo_cmd="sudo"
  DEBIAN_FRONTEND=noninteractive $sudo_cmd apt-get install -y texlive-xetex texlive-latex-recommended texlive-latex-extra texlive-fonts-recommended lmodern
fi
echo "Add to PATH: /opt/quarto/opt/quarto/bin:$ENV_DIR/bin"
