# Build command:
#   $ docker build -t aprokop/fortrilinos-stack:latest -f Dockerfile  .
ARG BASE=nvidia/cuda:11.0.3-devel-ubuntu18.04
FROM $BASE

ARG NPROC=8

RUN if test ${NV_CUDA_LIB_VERSION}; then apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/3bf863cc.pub; fi

# Avoid tzdata dialog
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -yq \
        autoconf \
        bc \
        build-essential \
        ccache \
        curl \
        environment-modules \
        gawk \
        gfortran \
        git \
        lcov \
        libatlas-base-dev \
        libbz2-dev \
        python2.7-dev \
        tmux \
        unzip \
        valgrind \
        vim \
        wget \
        && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install CMake
ENV CMAKE_DIR=/opt/cmake
RUN CMAKE_VERSION=3.16.4 && \
    CMAKE_URL=https://github.com/Kitware/CMake/releases/download/v${CMAKE_VERSION} && \
    CMAKE_SCRIPT=cmake-${CMAKE_VERSION}-Linux-x86_64.sh && \
    CMAKE_SHA256=cmake-${CMAKE_VERSION}-SHA-256.txt && \
    SCRATCH_DIR=/scratch && mkdir -p ${SCRATCH_DIR} && cd ${SCRATCH_DIR} && \
    wget --quiet ${CMAKE_URL}/${CMAKE_SHA256} && \
    wget --quiet ${CMAKE_URL}/${CMAKE_SHA256}.asc && \
    wget --quiet ${CMAKE_URL}/${CMAKE_SCRIPT} && \
    grep ${CMAKE_SCRIPT} ${CMAKE_SHA256} | sha256sum --check && \
    mkdir -p ${CMAKE_DIR} && \
    sh ${CMAKE_SCRIPT} --skip-license --prefix=${CMAKE_DIR} && \
    rm -rf ${SCRATCH_DIR}
ENV PATH=${CMAKE_DIR}/bin:$PATH

# Install OpenMPI
ENV OPENMPI_DIR=/opt/openmpi
RUN OPENMPI_VERSION=4.0.3 && \
    OPENMPI_VERSION_SHORT=$(echo "$OPENMPI_VERSION" | cut -d. -f1,2) && \
    OPENMPI_SHA1=d958454e32da2c86dd32b7d557cf9a401f0a08d3 && \
    OPENMPI_URL=https://download.open-mpi.org/release/open-mpi/v${OPENMPI_VERSION_SHORT}/openmpi-${OPENMPI_VERSION}.tar.bz2 && \
    OPENMPI_ARCHIVE=openmpi-${OPENMPI_VERSION}.tar.bz2 && \
    SCRATCH_DIR=/scratch && mkdir -p ${SCRATCH_DIR} && cd ${SCRATCH_DIR} && \
    wget --quiet ${OPENMPI_URL} --output-document=${OPENMPI_ARCHIVE} && \
    echo "${OPENMPI_SHA1} ${OPENMPI_ARCHIVE}" | sha1sum -c && \
    mkdir -p openmpi && \
    tar -xf ${OPENMPI_ARCHIVE} -C openmpi --strip-components=1 && \
    mkdir -p build && cd build && \
    ../openmpi/configure --prefix=${OPENMPI_DIR} ${CUDA_OPTIONS} CFLAGS=-w && \
    make -j${NPROCS} install && \
    rm -rf ${SCRATCH_DIR}
ENV PATH=${OPENMPI_DIR}/bin:$PATH

# Workaround for Kokkos to find libcudart
ENV LD_LIBRARY_PATH=/usr/local/cuda/targets/x86_64-linux/lib:${LD_LIBRARY_PATH}

# Install Trilinos (13.0.0)
ENV TRILINOS_DIR=/opt/trilinos
RUN export TRILINOS_HASH=9fec35276d846a667bc668ff4cbdfd8be0dfea08 && \
    export TRILINOS_URL=https://github.com/trilinos/Trilinos/archive/${TRILINOS_HASH}.tar.gz && \
    export TRILINOS_ARCHIVE=trilinos-${TRILINOS_HASH}.tar.gz && \
    SCRATCH_DIR=/scratch && mkdir -p ${SCRATCH_DIR} && cd ${SCRATCH_DIR} && \
    wget --quiet ${TRILINOS_URL} --output-document=${TRILINOS_ARCHIVE} && \
    mkdir -p trilinos && \
    tar -xf ${TRILINOS_ARCHIVE} -C trilinos --strip-components=1 && \
    mkdir -p build && cd build && \
    cmake \
        -D CMAKE_BUILD_TYPE=RelWithDebInfo \
        -D BUILD_SHARED_LIBS=ON \
        \
        -D CMAKE_INSTALL_PREFIX=${TRILINOS_DIR} \
        \
        -D TPL_ENABLE_Boost=OFF \
              -D TPL_ENABLE_CUDA=ON \
        -D TPL_ENABLE_CUDA=OFF \
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
    make -j${NPROC} install && \
    rm -rf ${SCRATCH_DIR}
ENV TRILINOS_DIR=/opt/trilinos
