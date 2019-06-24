# Build command:
#   $ docker build -t aprokop/fortrilinos-stack:latest -f Dockerfile_stack  .
FROM ubuntu:16.04

ARG NPROC=14

# 1. Need gcc-7 to install pgmath using Spack (required by Flang)
# 2. Do not install cmake as we'll need cmake 3.10+ for Trilinos
# 3. Need gawk for pgmath as default mawk segfaults
# 4. Need bc for kokkos build system (and other misc trilinos build)
RUN apt-get update && apt-get install -y \
        autoconf \
        bc \
        build-essential \
        curl \
        environment-modules \
        gfortran \
        git \
        lcov \
        libatlas-base-dev \
        libbz2-dev \
        libfreetype6-dev \
        libopenmpi-dev \
        libpng-dev \
        libsqlite3-dev \
        libssl-dev \
        libxft-dev \
        python2.7-dev \
        tmux \
        unzip \
        valgrind \
        vim \
        wget \
        zlib1g-dev \
        && \
    apt-get install -y gawk && \
    apt-get install -y software-properties-common && \
    add-apt-repository -y ppa:ubuntu-toolchain-r/test && \
    apt-get update && apt-get install -y \
        gcc-7 \
        g++-7 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* && \
    ln -s /usr/bin/python2.7 /usr/bin/python

ENV PREFIX=/scratch
RUN mkdir -p ${PREFIX} && \
    cd ${PREFIX} && \
    mkdir archive && \
    mkdir source && \
    mkdir build && \
    mkdir install

# append the option flag --allow-run-as-root to mpiexec
RUN echo '#!/usr/bin/env bash' > /usr/local/bin/mpiexec && \
    echo '/usr/bin/mpiexec --allow-run-as-root "$@"' >> /usr/local/bin/mpiexec && \
    chmod +x /usr/local/bin/mpiexec

ENV HOME=/root

# Install CMake
ENV CMAKE_DIR=/opt/cmake
RUN CMAKE_VERSION=3.13.4 && \
    CMAKE_VERSION_SHORT=3.13 && \
    CMAKE_URL=https://cmake.org/files/v${CMAKE_VERSION_SHORT}/cmake-${CMAKE_VERSION}-Linux-x86_64.sh && \
    CMAKE_SCRIPT=cmake-${CMAKE_VERSION}-Linux-x86_64.sh && \
    wget --quiet ${CMAKE_URL} --output-document=${CMAKE_SCRIPT} && \
    mkdir -p ${CMAKE_DIR} && \
    sh ${CMAKE_SCRIPT} --skip-license --prefix=${CMAKE_DIR} && \
    rm ${CMAKE_SCRIPT}
ENV PATH=${CMAKE_DIR}/bin:$PATH

# Install Flang
ENV FLANG_DIR=/opt/flang
RUN FLANG_VERSION=20190329 && \
    FLANG_ARCHIVE=flang-${FLANG_VERSION}-x86-70.tgz && \
    FLANG_URL=https://github.com/flang-compiler/flang/releases/download/flang_${FLANG_VERSION}/${FLANG_ARCHIVE} && \
    wget --quiet ${FLANG_URL} --output-document=${FLANG_ARCHIVE} && \
    mkdir -p ${FLANG_DIR} && \
    tar -xvf ${FLANG_ARCHIVE} -C ${FLANG_DIR} && \
    rm ${FLANG_ARCHIVE}
ENV PATH=${FLANG_DIR}/bin:$PATH
ENV LD_LIBRARY_PATH=${FLANG_DIR}/lib:$LD_LIBRARY_PATH

# Download Trilinos (version specified in TrilinosVersion.cmake)
COPY "trilinos_version" "${PREFIX}/"
RUN export TRILINOS_HASH="$(cat ${PREFIX}/trilinos_version)" && \
    export TRILINOS_VERSION="${TRILINOS_HASH}" && \
    export TRILINOS_URL="https://github.com/trilinos/Trilinos/archive/${TRILINOS_HASH}.tar.gz" && \
    export TRILINOS_ARCHIVE="${PREFIX}/archive/trilinos-${TRILINOS_HASH}.tar.gz" && \
    export TRILINOS_SOURCE_DIR="${PREFIX}/source/trilinos/${TRILINOS_HASH}" && \
    export TRILINOS_BUILD_DIR="${PREFIX}/build/trilinos/${TRILINOS_HASH}" && \
    wget --quiet "${TRILINOS_URL}" --output-document="${TRILINOS_ARCHIVE}" && \
    mkdir -p "${TRILINOS_SOURCE_DIR}" && \
    tar -xf "${TRILINOS_ARCHIVE}" -C "${TRILINOS_SOURCE_DIR}" --strip-components=1 && \
    ln -s "${TRILINOS_SOURCE_DIR}" "${PREFIX}/source/trilinos/release" && \
    mkdir -p "${TRILINOS_BUILD_DIR}" && \
    rm -rf "${TRILINOS_ARCHIVE}"

ENV TRILINOS_DIR="/scratch/source/trilinos/release"
