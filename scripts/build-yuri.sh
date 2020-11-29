#!/bin/sh -e
BUILDSCRIPT_DIR="$(cd "$(dirname $BASH_SOURCE[0])" && pwd)"
SOURCE_DIR="$(cd "${BUILDSCRIPT_DIR}" && git rev-parse --show-toplevel)"

printf "\e[2;37mBuilding from ${SOURCE_DIR}\e[0m\n"
cd $SOURCE_DIR
mkdir build 2>/dev/null || true
cd build

module purge
SPACK_VIEW=$SPACK_ROOT/var/spack/environments/fortrilinos/.spack-env/view
export CMAKE_PREFIX_PATH=$SPACK_VIEW:$CMAKE_PREFIX_PATH
export PATH=$SPACK_VIEW/bin:$PATH

cmake -G Ninja \
  -DCMAKE_INSTALL_PREFIX:PATH=$SOURCE_DIR/install \
  -DCMAKE_BUILD_TYPE:STRING=Debug \
  -DBUILD_SHARED_LIBS:BOOL=ON \
  -DForTrilinos_EXAMPLES:BOOL=ON \
  -DForTrilinos_TESTING:BOOL=ON \
  -DForTrilinos_USE_SWIG_Fortran:BOOL=ON \
  ..
ninja -v
ctest --output-on-failure
