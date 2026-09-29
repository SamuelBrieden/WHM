#!/usr/bin/env bash
#
# setup_whm_env.sh — reproducible setup of the `whm` conda environment for the
# Web-Halo Model (WHM) project.
#
# Builds NATIVE arm64 by default (miniforge / mambaforge on Apple Silicon).
# Falls back to x86_64 / Rosetta (CONDA_SUBDIR=osx-64) ONLY if the
# baccoemu/tensorflow stack fails to import on arm64.
#
# Verified working architecture (2026-06-02): arm64 (osx-arm64).
#   python 3.11.15, numpy 2.x, scipy 1.17.x, tensorflow 2.21.0, baccoemu 2.3.0,
#   camb 1.4.0 (LOCAL modified WHM-CAMB), euclidemu2, velocileptors.
#
# Usage:
#   bash setup_whm_env.sh            # native arm64 (default)
#   WHM_ARCH=osx-64 bash setup_whm_env.sh   # force x86_64 / Rosetta fallback
#
set -euo pipefail

ENV_NAME="${ENV_NAME:-whm}"
PY_VER="${PY_VER:-3.11}"
WHM_ARCH="${WHM_ARCH:-osx-arm64}"   # osx-arm64 (default) | osx-64 (fallback)

# Resolve repo root: this script lives in WHM-notebooks/env/ ; WHM-CAMB is the sibling folder ../../WHM-CAMB
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CAMB_DIR="${WHM_CAMB_DIR:-$(cd "${SCRIPT_DIR}/../.." && pwd)/WHM-CAMB}"

# Prefer a native-arm64 miniforge/mambaforge install for the `conda` binary.
CONDA_BIN="$(command -v conda || true)"
if [ -x "$HOME/miniforge3/bin/conda" ]; then
  CONDA_BIN="$HOME/miniforge3/bin/conda"
elif [ -x "$HOME/mambaforge/bin/conda" ]; then
  CONDA_BIN="$HOME/mambaforge/bin/conda"
fi
if [ -z "${CONDA_BIN}" ]; then
  echo "ERROR: no conda found. Install miniforge (arm64) first:" >&2
  echo "  https://github.com/conda-forge/miniforge#install" >&2
  exit 1
fi
CONDA_ROOT="$("${CONDA_BIN}" info --base)"
# shellcheck disable=SC1091
source "${CONDA_ROOT}/etc/profile.d/conda.sh"

echo "==> conda: ${CONDA_BIN}  (base ${CONDA_ROOT})"
echo "==> target arch: ${WHM_ARCH}   env: ${ENV_NAME}   python ${PY_VER}"

# ---------------------------------------------------------------------------
# 1. Create the conda env from conda-forge (override defaults channel).
# ---------------------------------------------------------------------------
CONDA_PKGS=(
  "python=${PY_VER}" numpy scipy sympy matplotlib jupyter cython
  pyfftw fftw gsl gfortran make cmake pkg-config c-compiler cxx-compiler
)

if [ "${WHM_ARCH}" = "osx-64" ]; then
  echo "==> Creating env as x86_64 (Rosetta fallback)."
  CONDA_SUBDIR=osx-64 "${CONDA_BIN}" create -y -n "${ENV_NAME}" \
    --override-channels -c conda-forge "${CONDA_PKGS[@]}"
  conda activate "${ENV_NAME}"
  conda config --env --set subdir osx-64
else
  echo "==> Creating env as native arm64."
  "${CONDA_BIN}" create -y -n "${ENV_NAME}" \
    --override-channels -c conda-forge "${CONDA_PKGS[@]}"
  conda activate "${ENV_NAME}"
fi

# setuptools 81+ removed pkg_resources, which WHM-CAMB/camb/_compilers.py imports.
# Pin <81 so the (no-isolation) build of the local CAMB succeeds.
pip install "setuptools<81" wheel

# ---------------------------------------------------------------------------
# 2. GATE: TensorFlow + baccoemu. This is the only real architecture risk.
#    Run it FIRST, before building CAMB/CLASS.
# ---------------------------------------------------------------------------
echo "==> GATE: installing tensorflow + baccoemu and instantiating Matter_powerspectrum()"
pip install tensorflow baccoemu
if ! python -c "import baccoemu; baccoemu.Matter_powerspectrum(); print('bacco OK')"; then
  echo "" >&2
  echo "ERROR: baccoemu/tensorflow failed to import on ${WHM_ARCH}." >&2
  if [ "${WHM_ARCH}" != "osx-64" ]; then
    echo "Retry under Rosetta with:" >&2
    echo "  conda env remove -n ${ENV_NAME} -y && WHM_ARCH=osx-64 bash ${BASH_SOURCE[0]}" >&2
  fi
  exit 1
fi

# ---------------------------------------------------------------------------
# 3. Build the LOCAL modified WHM-CAMB from source (compiles fortran/ via
#    gfortran). NEVER use pip's published camb. --no-build-isolation so the
#    pinned setuptools<81 (with pkg_resources) is used by setup.py.
# ---------------------------------------------------------------------------
echo "==> Building local WHM-CAMB at ${CAMB_DIR}"
rm -f "${CAMB_DIR}/camb/camblib.so" "${CAMB_DIR}/build/lib/camb/camblib.so" 2>/dev/null || true
make -C "${CAMB_DIR}/fortran" clean 2>/dev/null || true
pip install -e "${CAMB_DIR}" --no-build-isolation
python -c "import camb; from camb.nonlinear import halofit_brieden2025_tweaked; \
print('camb', camb.__version__, 'brieden modes OK', camb.__file__)"

# ---------------------------------------------------------------------------
# 4. Remaining deps.
# ---------------------------------------------------------------------------
echo "==> Installing velocileptors, euclidemu2, FLAMINGO baryon response emulator"
pip install git+https://github.com/sfschen/velocileptors
pip install euclidemu2
pip install git+https://github.com/FLAMINGOSIM/FlamingoBaryonResponseEmulator

# ---------------------------------------------------------------------------
# 5. CLASS (classy) is built in a SEPARATE hand-off (class_public). Not here.
# ---------------------------------------------------------------------------
echo "==> NOTE: CLASS / classy is NOT built by this script (separate hand-off)."

# ---------------------------------------------------------------------------
# 6. Smoke-test.
# ---------------------------------------------------------------------------
echo "==> Running import smoke-test"
python "${SCRIPT_DIR}/import_smoketest.py"

echo "==> Done. Activate with: conda activate ${ENV_NAME}"
