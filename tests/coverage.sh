#!/usr/bin/env bash
#
# coverage.sh -- build APFEL++ with coverage instrumentation, run the full
# regression suite through ctest, and produce a code-coverage report.
#
# Usage (from anywhere):
#     tests/coverage.sh
#
# The report is written to:
#     build-coverage/coverage-report/index.html   (clickable, per-file)
#     build-coverage/coverage.xml                  (Cobertura, for CI)
# and a summary is printed to the terminal.
#
# Requirements: cmake, a C++ compiler, and gcovr (`pip install gcovr`).
# On macOS the Apple-Clang gcov data is parsed via `xcrun llvm-cov gcov`;
# on Linux the system `gcov` is used. Both are selected automatically.

set -euo pipefail

# Resolve repository root (this script lives in <root>/tests).
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
BUILD_DIR="${ROOT_DIR}/build-coverage"
REPORT_DIR="${BUILD_DIR}/coverage-report"

# --- Sanity check: gcovr must be installed -------------------------------
if ! command -v gcovr >/dev/null 2>&1; then
    cat >&2 <<'EOF'
ERROR: gcovr was not found on PATH.

Install it, e.g. in a virtual environment:
    python3 -m venv ~/.venvs/gcovr
    ~/.venvs/gcovr/bin/pip install gcovr
    export PATH="$HOME/.venvs/gcovr/bin:$PATH"

then re-run tests/coverage.sh
EOF
    exit 1
fi

# --- Select the gcov backend per platform --------------------------------
GCOV_ARGS=()
if [[ "$(uname -s)" == "Darwin" ]]; then
    # Apple's /usr/bin/gcov cannot read Clang coverage data; use llvm-cov.
    GCOV_ARGS+=(--gcov-executable "xcrun llvm-cov gcov")
fi

JOBS="$( (command -v nproc >/dev/null && nproc) || sysctl -n hw.ncpu || echo 4)"

echo "==> Configuring instrumented build in ${BUILD_DIR}"
cmake -S "${ROOT_DIR}" -B "${BUILD_DIR}" -DAPFEL_COVERAGE=ON

echo "==> Building (-j${JOBS})"
cmake --build "${BUILD_DIR}" -j"${JOBS}"

echo "==> Running regression suite via ctest"
# Keep going even if a test fails: a failing test still produces coverage data,
# and we want the report regardless. ctest's own status is reported at the end.
ctest_status=0
( cd "${BUILD_DIR}" && ctest -j"${JOBS}" --output-on-failure ) || ctest_status=$?

echo "==> Generating coverage report with gcovr"
mkdir -p "${REPORT_DIR}"
gcovr \
    --root "${ROOT_DIR}" \
    --object-directory "${BUILD_DIR}" \
    "${GCOV_ARGS[@]}" \
    --gcov-ignore-errors=all \
    --gcov-ignore-parse-errors=all \
    --filter "${ROOT_DIR}/src/" \
    --filter "${ROOT_DIR}/inc/" \
    --exclude-unreachable-branches \
    --html-details "${REPORT_DIR}/index.html" \
    --xml "${BUILD_DIR}/coverage.xml" \
    --print-summary

echo
echo "==> HTML report: ${REPORT_DIR}/index.html"
echo "==> Cobertura  : ${BUILD_DIR}/coverage.xml"

if [[ "${ctest_status}" -ne 0 ]]; then
    echo "WARNING: ctest reported failures (exit ${ctest_status}); coverage above still reflects what ran." >&2
fi
exit "${ctest_status}"
