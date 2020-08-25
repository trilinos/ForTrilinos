# Build command:
#   $ docker build -t sethrj/fortrilinos-stack:latest -f Dockerfile  .
FROM aprokop/fortrilinos-stack:latest

# Install Trilinos
RUN cd /scratch/build/trilinos \
    && cmake \
      -DTrilinos_ENABLE_Amesos2:BOOL=ON \
      -DTrilinos_ENABLE_Anasazi:BOOL=ON \
      -DTrilinos_ENABLE_Belos:BOOL=ON \
      -DTrilinos_ENABLE_Epetra:BOOL=OFF \
      -DTrilinos_ENABLE_Ifpack2:BOOL=ON \
      -DTrilinos_ENABLE_MueLu:BOOL=ON \
      -DTrilinos_ENABLE_NOX:BOOL=ON \
      -DTrilinos_ENABLE_Stratimikos:BOOL=ON \
      -DTrilinos_ENABLE_Teuchos:BOOL=ON \
      -DTrilinos_ENABLE_Thyra:BOOL=ON \
      -DTrilinos_ENABLE_Tpetra:BOOL=ON \
      -DTpetra_INST_INT_INT:BOOL=OFF \
      -DTpetra_INST_INT_LONG:BOOL=OFF \
      -DTpetra_INST_INT_LONG_LONG:BOOL=ON \
      -DTpetra_INST_FLOAT:BOOL=OFF \
      -DTpetra_ENABLE_DEPRECATED_CODE=OFF \
      -DTPL_ENABLE_MPI:BOOL=ON \
      -DTPL_ENABLE_Boost:BOOL=OFF \
      -DTPL_ENABLE_CUDA:BOOL=OFF \
      -DTPL_ENABLE_DLlib:BOOL=OFF \
      -DTrilinos_ENABLE_DEBUG:BOOL=OFF \
      -DTrilinos_ENABLE_Fortran:BOOL=ON \
      -DTrilinos_ENABLE_EXAMPLES:BOOL=OFF \
      -DTrilinos_ENABLE_TESTS:BOOL=OFF \
      -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=OFF \
      -DTrilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES:BOOL=OFF \
      -DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF \
      -DTrilinos_ENABLE_SECONDARY_TESTED_CODE:BOOL=ON \
      -DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON \
      -DCMAKE_BUILD_TYPE:STRING=RelWithDebInfo \
      -DCMAKE_INSTALL_PREFIX:STRING=/scratch/install/trilinos \
      -DBUILD_SHARED_LIBS:BOOL=ON \
    ${TRILINOS_DIR} \
    && make -j${NPROC} install \
    && rm -rf /scratch/build/trilinos

