#!/usr/bin/env bash

# number of processes with default value
: ${NPROC:=8}
# cleanup workspace
cd ${TRILINOS_DIR}/packages/ForTrilinos
[ -d build ] && rm -rf build
mkdir build && cd build
# configure trilinos with fortrilinos
../scripts/docker_cmake -D Trilinos_ENABLE_COVERAGE_TESTING=ON
# build
# build twice, as the first one fails with
#       use TEST_CALLS_FILE
#           1
#   Fatal Error: Can't open module file epetra_blockmap_test_calls.mod for reading at (1): No such file or directory
make -j${NPROC} -i
make -j${NPROC} -i
# run the unit tests
ctest -j${NPROC} --no-compress-output -T Test
# collect coverage data
lcov --capture --directory packages/ForTrilinos --output-file lcov.info
# upload it to codecov
curl -s https://codecov.io/bash -o codecov_bash_uploader
chmod +x codecov_bash_uploader
./codecov_bash_uploader -Z -X gcov -f lcov.info
