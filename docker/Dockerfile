# Build command:
#   $ docker build -t sethrj/fortrilinos-stack:latest -f Dockerfile  .
ARG BASE=nvidia/cuda:11.2.2-devel-ubuntu20.04
FROM $BASE

ARG NPROC=8

# Avoid tzdata dialog
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get -yqq update && apt-get install -yq \
        autoconf \
        bc \
        build-essential \
        ca-certificates \
        ccache \
        curl \
        gawk \
        gfortran \
        g++ \
        git \
        lcov \
        libatlas-base-dev \
        libbz2-dev \
        libssl-dev \
        mpich \
        ninja-build \
        python3-dev \
        valgrind \
        vim \
        wget \
        && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

ENV CMAKE_DIR=/opt/cmake
RUN CMAKE_URL=https://github.com/Kitware/CMake/releases/download/v3.23.3/cmake-3.23.3.tar.gz && \
    CMAKE_ARCHIVE=cmake.tar.gz && \
    SCRATCH_DIR=/scratch && mkdir -p ${SCRATCH_DIR} && cd ${SCRATCH_DIR} && \
    wget --quiet ${CMAKE_URL}  --output-document=${CMAKE_ARCHIVE} && \
    mkdir -p ${CMAKE_DIR} && \
    tar -xf ${CMAKE_ARCHIVE} --strip-components=1 && \
    ./bootstrap --prefix=${CMAKE_DIR} --generator=Ninja && \
    ninja install && \
    rm -rf ${SCRATCH_DIR}
ENV PATH=${CMAKE_DIR}/bin:$PATH

# Install Trilinos
ENV TRILINOS_DIR=/opt/trilinos
RUN TRILINOS_VERSION=13-4-1 && \
    TRILINOS_URL=https://github.com/trilinos/Trilinos/archive/refs/tags/trilinos-release-${TRILINOS_VERSION}.tar.gz && \
    TRILINOS_ARCHIVE=trilinos-${TRILINOS_VERSION}.tar.gz && \
    SCRATCH_DIR=/scratch && mkdir -p ${SCRATCH_DIR} && cd ${SCRATCH_DIR} && \
    wget --quiet ${TRILINOS_URL} --output-document=${TRILINOS_ARCHIVE} && \
    mkdir -p trilinos && \
    tar -xf ${TRILINOS_ARCHIVE} -C trilinos --strip-components=1 && \
    mkdir -p build && cd build && \
    cmake \
        -G Ninja \
        -D CMAKE_BUILD_TYPE=RelWithDebInfo \
        -D BUILD_SHARED_LIBS=ON \
        \
        -D CMAKE_INSTALL_PREFIX=${TRILINOS_DIR} \
        \
        -D TPL_ENABLE_Boost=OFF \
        -D TPL_ENABLE_DLlib=OFF \
        -D TPL_ENABLE_MPI=ON \
        \
        -D Trilinos_ENABLE_Fortran=ON \
        \
        -D Trilinos_ENABLE_EXPLICIT_INSTANTIATION=ON \
        -D Trilinos_ENABLE_DEBUG=OFF \
        -D Trilinos_ENABLE_EXAMPLES=OFF \
        -D Trilinos_ENABLE_TESTS=OFF \
        -D Trilinos_ENABLE_ALL_PACKAGES=OFF \
        -D Trilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES=OFF \
        -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES=OFF \
        -D Trilinos_ENABLE_SECONDARY_TESTED_CODE=ON \
        \
        -D Trilinos_ENABLE_Amesos2=ON \
        -D Trilinos_ENABLE_Anasazi=ON \
        -D Trilinos_ENABLE_Belos=ON \
        -D Trilinos_ENABLE_Epetra=OFF \
        -D Trilinos_ENABLE_Ifpack2=ON \
        -D Trilinos_ENABLE_Kokkos=ON \
            -D Kokkos_ENABLE_Serial=ON \
            -D Kokkos_ENABLE_OpenMP=OFF \
            -D Kokkos_ENABLE_Cuda=OFF \
            -D Kokkos_ENABLE_Cuda_UVM=ON \
            -D Kokkos_ENABLE_Cuda_Lambda=ON \
            -D Kokkos_ARCH_VOLTA70=ON \
        -D Trilinos_ENABLE_MueLu=ON \
        -D Trilinos_ENABLE_NOX=ON \
        -D Trilinos_ENABLE_Stratimikos=ON \
        -D Trilinos_ENABLE_Teuchos=ON \
        -D Trilinos_ENABLE_Thyra=ON \
        -D Trilinos_ENABLE_Tpetra=ON \
            -D Tpetra_INST_SERIAL=ON \
            -D Tpetra_INST_INT_INT=OFF \
            -D Tpetra_INST_INT_LONG=OFF \
            -D Tpetra_INST_INT_LONG_LONG=ON \
            -D Tpetra_INST_FLOAT=OFF \
            -D Tpetra_ENABLE_DEPRECATED_CODE=OFF \
        \
    ../trilinos && \
    ninja -j${NPROC} install && \
    rm -rf ${SCRATCH_DIR}
ENV TRILINOS_DIR=/opt/trilinos
