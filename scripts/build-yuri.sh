#!/bin/sh -e
BUILDSCRIPT_DIR="$(cd "$(dirname $BASH_SOURCE[0])" && pwd)"
SOURCE_DIR="$(cd "${BUILDSCRIPT_DIR}" && git rev-parse --show-toplevel)"

printf "\e[2;37mBuilding from ${SOURCE_DIR}\e[0m\n"
cd $SOURCE_DIR
mkdir build 2>/dev/null || true
cd build

module purge
SPACK_VIEW=/usr/local/spack/var/spack/environments/fortrilinos/.spack-env/view
export CMAKE_PREFIX_PATH=$SPACK_VIEW:$CMAKE_PREFIX_PATH
export PATH=$SPACK_VIEW/bin:$PATH

cmake -G Ninja \
  -DCMAKE_INSTALL_PREFIX:PATH=$SOURCE_DIR/install \
  -DCMAKE_BUILD_TYPE:STRING=Debug \
  -DBUILD_SHARED_LIBS:BOOL=ON \
  -DForTrilinos_USE_SWIG_Fortran:BOOL=ON \
  -DSWIG_EXECUTABLE:FILENAME=/rnsdhpc/code/swig-old/swig \
  -DSWIG_DIR:FILENAME=/rnsdhpc/code/swig-old/Lib \
  ..
ninja -v
ctest --output-on-failure
