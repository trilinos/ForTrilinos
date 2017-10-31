#!/usr/bin/env bash

set -e

# number of processes with default value
: ${NPROC:=8}
# bind mount ForTrilinos source dir into Trilinos base dir
mkdir ${TRILINOS_DIR}/packages/ForTrilinos
mount --bind ${FORTRILINOS_DIR} ${TRILINOS_DIR}/packages/ForTrilinos
# cleanup workspace
cd ${TRILINOS_DIR}/packages/ForTrilinos
[ -d build ] && rm -rf build
mkdir build && cd build
# NOTE: relative paths are invalid after configuration when ForTrilinos source dir is
# not directly mounted into Trilinos base source dir. We build elsewhere and
# move the build directory afterwards...
# configure trilinos with fortrilinos

if [ "${BUILD_TYPE}" == "gcc54-mpi" ]; then
  ../scripts/docker_cmake -D Trilinos_ENABLE_COVERAGE_TESTING=ON
elif [ "${BUILD_TYPE}" == "gcc54-serial" ]; then
  ../scripts/docker_cmake_serial
else
    echo "Unknown BUILD_TYPE"
    exit 1
fi

# build
# build twice, as the first one fails with
#       use TEST_CALLS_FILE
#           1
#   Fatal Error: Can't open module file epetra_blockmap_test_calls.mod for reading at (1): No such file or directory
make -j${NPROC} -i
make -j${NPROC} -i
# run the unit tests
ctest -j${NPROC} --no-compress-output -T Test
# upload code coverage only once
if [ "${BUILD_TYPE}" == "gcc54-mpi"  ]
then
  # collect coverage data
  lcov --capture --directory packages/ForTrilinos --output-file lcov.info
  # upload it to codecov
  curl -s https://codecov.io/bash -o codecov_bash_uploader
  chmod +x codecov_bash_uploader
  ./codecov_bash_uploader -Z -X gcov -f lcov.info
fi
