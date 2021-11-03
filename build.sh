#!/usr/bin/env bash
cmake                               \
  -DCMAKE_BUILD_TYPE=Release        \
  -DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -O3 -g -march=znver2"  \
  -DINCLUDE_DIR=${1}  \
  -DLINK_DIR=${2}  \
..
