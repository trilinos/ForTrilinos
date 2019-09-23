#!/usr/bin/env bash

set -e

# number of processes with default value
: "${NPROC:=8}"
# get most recent version of develop
if [ -n "${TRILINOS_VERSION}" ]; then
  # Remove pre-existing Trilinos version in the container
  rm -rf "${TRILINOS_DIR}"
  # Download speficifed Trilinos version
  PREFIX="/scratch"
  TRILINOS_URL="https://github.com/trilinos/Trilinos/archive/${TRILINOS_VERSION}.tar.gz"
  TRILINOS_ARCHIVE="${PREFIX}/archive/trilinos-${TRILINOS_VERSION}.tar.gz"
  wget --quiet "${TRILINOS_URL}" --output-document="${TRILINOS_ARCHIVE}"
  mkdir "${TRILINOS_DIR}"
  tar -xf "${TRILINOS_ARCHIVE}" -C "${TRILINOS_DIR}" --strip-components=1
  rm -rf "${TRILINOS_ARCHIVE}"
  echo "Updating container version of Trilinos to \"${TRILINOS_VERSION}\""
else
  echo "Using container version of Trilinos"
fi
# bind mount ForTrilinos source dir into Trilinos base dir
mkdir "${TRILINOS_DIR}/ForTrilinos"
mount --bind "${FORTRILINOS_DIR}" "${TRILINOS_DIR}/ForTrilinos"
# cleanup workspace
cd "${TRILINOS_DIR}/ForTrilinos"
[ -d build ] && rm -rf build
mkdir build && cd build
# NOTE: relative paths are invalid after configuration when ForTrilinos source dir is
# not directly mounted into Trilinos base source dir. We build elsewhere and
# move the build directory afterwards...
# configure trilinos with fortrilinos

if [ "${BUILD_TYPE}" == "gcc74-mpi" ]; then
  ../scripts/docker_cmake -D Trilinos_ENABLE_COVERAGE_TESTING=ON

elif [ "${BUILD_TYPE}" == "gcc74-mpi-openmp" ]; then
  ../scripts/docker_cmake \
    -D Trilinos_ENABLE_OpenMP=ON \
    -D Kokkos_ENABLE_OpenMP=ON \
    -D Tpetra_INST_OPENMP=ON

elif [ "${BUILD_TYPE}" == "gcc74-mpi-cuda" ]; then
  export OMPI_CXX="${TRILINOS_DIR}/packages/kokkos/bin/nvcc_wrapper"
  ../scripts/docker_cmake \
    -D TPL_ENABLE_CUDA=ON \
    -D Trilinos_CXX11_FLAGS="-std=c++11 -expt-extended-lambda" \
    -D KOKKOS_ARCH="Volta70" \
    -D Kokkos_ENABLE_Cuda=ON \
    -D Kokkos_ENABLE_Cuda_UVM=ON \
    -D Kokkos_ENABLE_Cuda_Lambda=ON \
    -D Tpetra_INST_CUDA=ON

elif [ "${BUILD_TYPE}" == "flang70-serial" ]; then
  source ../scripts/docker_flang70_env.sh
  # - No MPI with Flang
  #   The system installation of openmpi produces incompatible .mod files.
  # - Use "RelWithDebInfo" build type
  #   Flang does not compile in Debug mode. It produces spurious undefined
  #   symbols ending at _tbp_, resulting in undefined references during
  #   linking stage. This does not happen in RelWithDebInfo mode.
  ../scripts/docker_cmake_serial -DCMAKE_BUILD_TYPE="RelWithDebInfo"

else
  echo "Unknown BUILD_TYPE"
  exit 1
fi

# build
make -j"${NPROC}" -i
# run the unit tests
ctest -j"${NPROC}" --no-compress-output --output-on-failure -T Test
# upload code coverage only once
if [ "${BUILD_TYPE}" == "gcc74-mpi"  ]; then
  # collect coverage data
  lcov --capture --directory ForTrilinos --output-file lcov.info
  # upload it to codecov
  curl -s https://codecov.io/bash -o codecov_bash_uploader
  chmod +x codecov_bash_uploader
  ./codecov_bash_uploader -Z -X gcov -f lcov.info
fi
