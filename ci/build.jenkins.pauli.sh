#!/usr/bin/env bash

# This script expects the following environment variables
# COMPILER - the compiler to build with
# BUILD_CONFIG - 'Debug' | 'Release' | 'RelWithDebInfo'
# BUILD_DIR - (optional) A directory name to build in (default: formed from COMPILER and MPI_VERSION)
# PHASE - 'cmake' | 'make' | 'test' (all of the above if empty)
# The scripts creates directories under ./build.tmp/

BASE_DIR=$(/bin/pwd)

# This function sets build environment unless it's already set (as determined by env. var `build_environment_set`)
function setup_environment() {
  [[ $build_environment_set == 1 ]] && return 0

  # Jenkins doesn't seem to load environment by default
  . /etc/profile

  module purge

  # Replace '_' by '/'
  _COMPILER_MODULE="BuildEnv/${COMPILER/_/-}"
  _COMPILER_PREFIX="${COMPILER/_//}"

  case $COMPILER in
  gcc_10.2.0)
    module add ${_COMPILER_MODULE}
    ;;
  llvm_5.0.1)
    module add ${_COMPILER_MODULE}
    ;;
  intel_19.0.2.187)
    module add ${_COMPILER_MODULE}
    ;;
  *)
    echo "Unsupported compiler passed via COMPILER='$COMPILER'; valid values are:" 2>&1
    echo "gcc_7.3.0 llvm_5.0.1 intel_19.0.2.187"
    exit 1
    ;;

  esac

  case $BUILD_CONFIG in
  Debug)
    echo "Performing Debug build"
    ;;
  Release)
    echo "Performing Release build"
    ;;
  RelWithDebInfo)
    echo "Performing Release build with debug info"
    ;;
  *)
    echo "Unsupported build type passed via BUILD_CONFIG='$BUILD_CONFIG'; valid values are:" 2>&1
    echo "Debug Release RelWithDebInfo"
    exit 1
    ;;
  esac

  module add openmpi/${_COMPILER_PREFIX}/3.1.4
  module add alpscore/${_COMPILER_PREFIX}/2.3.1a-mpi

  [[ $BUILD_DIR ]] || BUILD_DIR="build.tmp/${COMPILER}"

  mkdir -pv "$BUILD_DIR"
  cd "$BUILD_DIR"

  build_environment_set=1
}

function run_cmake() {
  rm -rf *
  cmake "$BASE_DIR" \
    -DCMAKE_BUILD_TYPE=$BUILD_CONFIG \
    -DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -O3 -g -march=znver2"  \
    -DINCLUDE_DIR=/data/tblommel/Install/include  \
    -DLINK_DIR=/data/tblommel/Install/lib  \
  ..
}

function run_make() {
  make -j8
}

function run_test() {
  make test
}

build_environment_set=0

if [[ $PHASE == 'cmake' || $PHASE == '' ]]; then
  setup_environment
  run_cmake || exit $?
fi

if [[ $PHASE == 'make' || $PHASE == '' ]]; then
  setup_environment
  run_make || exit $?
fi

if [[ $PHASE == 'test' || $PHASE == '' ]]; then
  setup_environment
  run_test || exit $?
fi

