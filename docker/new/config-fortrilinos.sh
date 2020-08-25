#!/bin/sh

SOURCE_DIR=/fortrilinos/src
BUILD_DIR=/scratch/build/fortrilinos
INSTALL_DIR=/scratch/install/fortrilinos
mkdir $BUILD_DIR
cd ${BUILD_DIR}

export CMAKE_INSTALL_PREFIX=/scratch/install/trilinos:$CMAKE_INSTALL_PREFIX

# Note that build type can be different from the trilinos build! That means
# debuggable wrappers but with optimized solvers.
cmake \
  -DCMAKE_BUILD_TYPE:STRING=Debug \
  -DForTrilinos_TESTING:BOOL=ON \
  -DCMAKE_INSTALL_PREFIX:STRING=$INSTALL_DIR \
  $SOURCE_DIR
make -j6
ctest -j6 --output-on-failure
make install
