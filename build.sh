#!/usr/bin/env bash
cmake                               \
  -DCMAKE_BUILD_TYPE=Release        \
  -DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -O3 -std=c++11"  \
  -DINCLUDE_DIR=~/Libraries/Install/include/  \
  -DLINK_DIR=~/Libraries/Install/lib  \
  -DDATA_DIR=/pauli-storage/tblommel/NEdyson_data/ \
  -DTIMING_DATA_DIR=/pauli-storage/tblommel/NEdyson_data/timing_data/ \
..
