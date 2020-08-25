#!/bin/sh -ex
###############################################################################
# File  : docker/test-installation.sh
###############################################################################

test -n "${FORTRILINOS_INSTALL_DIR}"

FORTRILINOS_INSTALL_DIR="$(cd ${FORTRILINOS_INSTALL_DIR} && pwd)"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

cd "${SCRIPT_DIR}/test-installation"
git clean -fxd .
mkdir build
cd build
export CMAKE_PREFIX_PATH=${FORTRILINOS_INSTALL_DIR}:${CMAKE_PREFIX_PATH}
cmake -G Ninja \
  -D CMAKE_INSTALL_PREFIX=${SCRIPT_DIR}/install \
  ..
ninja
ctest --output-on-failure

###############################################################################
# end of ci/scripts/test-installation.sh
###############################################################################
