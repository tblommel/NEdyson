#!/usr/bin/env bash
cmake                               \
  -DCMAKE_BUILD_TYPE=Release        \
  -DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -O3 -std=c++11"  \
  -DINCLUDE_DIR=/home/tblommel/Libraries/Install/include/  \
<<<<<<< HEAD
  -DLINK_DIR=/home/tblommel/Libraries/Install/lib  \
=======
  -DLINK_DIR=/home/tblommel/Libraries/Install/lib/  \
>>>>>>> 7210efc6c26be9cfc0cdcf199a20ae3fc06e7e93
..
