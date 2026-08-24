#!/usr/bin/env bash

set -euo pipefail

repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
tests_dir="$repo_root/tests"
htslib_dir="$repo_root/inst/include/htslib-1.18"
htslib_archive="$htslib_dir/libhts.a"
jobs="${JOBS:-4}"
cuda="${CUDA:-0}"
extra_libs=""

if [[ "$cuda" != "0" && "$cuda" != "1" ]]; then
    echo "CUDA must be 0 or 1; got: $cuda" >&2
    exit 2
fi
if [[ "$cuda" == "1" ]]; then
    extra_libs=" -lcudart"
fi

if [[ ! "$jobs" =~ ^[1-9][0-9]*$ ]]; then
    echo "JOBS must be a positive integer; got: $jobs" >&2
    exit 2
fi

if [[ ! -f "$htslib_archive" ]]; then
    echo "Building bundled HTSlib static archive..."
    if [[ ! -f "$htslib_dir/config.mk" ]]; then
        (
            cd "$htslib_dir"
            ./configure --disable-libcurl --without-libdeflate
        )
    fi
    make -C "$htslib_dir" -j"$jobs" libhts.a
fi

echo "Building and running the C++ test suite with $jobs job(s)..."
make -C "$tests_dir" -j"$jobs" -B test \
    CUDA="$cuda" \
    LDFLAGS= \
    "LIBS=$htslib_archive -llzma -lbz2 -lm -lz -lpthread$extra_libs"
