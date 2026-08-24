#!/usr/bin/env bash
# Antibody humanization pipeline - server installation
#
#   bash install.sh [--full] [--env-name humanize]
#
# Default (portable core): Python + standard library only, germline data
# bundled. Works on any Linux/macOS machine.
#   --full  additionally installs ANARCI (exact numbering cross-validation)
#           and IgBLAST binaries + NCBI germline database.
#
# ProteinMPNN and AlphaFold3 are NOT installed here: they are heavy
# GPU tools, usually deployed separately. Point the pipeline at them with
#   --mpnn-script /path/protein_mpnn.py --af3-mode local --af3-binary ...
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_NAME="humanize"
FULL=0

for arg in "$@"; do
  case "$arg" in
    --full) FULL=1 ;;
    --env-name) ENV_NAME="${2:-humanize}" ;;
  esac
done

PYTHON="${PYTHON:-python3}"

echo "==> Installing humanization pipeline into conda env '${ENV_NAME}'"

if command -v conda >/dev/null 2>&1; then
  conda create -n "${ENV_NAME}" python=3.11 -y
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate "${ENV_NAME}"
  PYTHON="$(which python)"
elif command -v mamba >/dev/null 2>&1; then
  mamba create -n "${ENV_NAME}" python=3.11 -y
  source "$(mamba info --base)/etc/profile.d/conda.sh"
  mamba activate "${ENV_NAME}"
  PYTHON="$(which python)"
else
  echo "    (conda not found; using ${PYTHON})"
fi

# --- optional exact-numbering cross-validation ---------------------------------
if [ "${FULL}" = "1" ]; then
  if command -v conda >/dev/null 2>&1; then
    conda install -y -n "${ENV_NAME}" -c bioconda -c conda-forge hmmer==3.3.2
  fi
  "${PYTHON}" -m pip install anarci abnumber biopython
else
  "${PYTHON}" -m pip install biopython 2>/dev/null || true
fi

# --- IgBLAST + NCBI germline (full mode) ---------------------------------------
if [ "${FULL}" = "1" ]; then
  echo "==> Downloading NCBI IgBLAST (binary + human germline FASTA)"
  "${PYTHON}" "${ROOT}/scripts/humanize/cli.py" setup-germline \
    --dir "${ROOT}/data/germline"
  echo "==> IgBLAST done (bundled germline remains available as fallback)"
else
  echo "==> Using bundled human germline data (373 V + 28 J genes, Kabat space)"
fi

echo
echo "==> Verify:"
echo "    ${PYTHON} ${ROOT}/scripts/humanize/cli.py setup-check"
echo "    ${PYTHON} ${ROOT}/tests/test_pipeline.py"
echo
echo "==> Run:"
echo "    ${PYTHON} ${ROOT}/scripts/humanize/cli.py run --input my_antibody.fasta --outdir outputs"
