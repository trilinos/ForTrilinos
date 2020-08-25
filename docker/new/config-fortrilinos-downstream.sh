#!/bin/sh

SOURCE_DIR=/fortrilinos/src/example/test-installation
BUILD_DIR=/scratch/build/fortrilinos-test
INSTALL_DIR=/scratch/install/fortrilinos-test
mkdir $BUILD_DIR
cd ${BUILD_DIR}

export CMAKE_INSTALL_PREFIX=/scratch/install/fortrilinos:$CMAKE_INSTALL_PREFIX

cmake \
  -DCMAKE_BUILD_TYPE:STRING=Debug \
  -DCMAKE_INSTALL_PREFIX:STRING=$INSTALL_DIR \
  $SOURCE_DIR
make -j6
ctest -j6 --output-on-failure
make install
