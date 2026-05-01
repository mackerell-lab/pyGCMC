#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${GCMC_BUILD_DIR:-$ROOT/build}"
BUILD_TYPE="${CMAKE_BUILD_TYPE:-Release}"
JOBS="${CMAKE_BUILD_PARALLEL_LEVEL:-4}"

cmake -S "$ROOT" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE="$BUILD_TYPE"
cmake --build "$BUILD_DIR" --target gcmc_cpu -j "$JOBS"

echo "$BUILD_DIR/bin/gcmc_cpu"
