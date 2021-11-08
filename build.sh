#!/usr/bin/env bash
cmake                               \
  -DCMAKE_BUILD_TYPE=Debug        \
  -DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -O3 -g -march=znver2"  \
  -DINCLUDE_DIR=~/Libraries/Install/include/  \
  -DLINK_DIR=~/Libraries/Install/lib  \
..
