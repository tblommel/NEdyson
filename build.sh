#!/usr/bin/env bash
cmake                               \
  -DCMAKE_BUILD_TYPE=Release        \
  -DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -O3 -std=c++11"  \
  -DINCLUDE_DIR=/usr/local/include/  \
  -DLINK_DIR=/usr/local/lib  \
..
