INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/packages/epetra/src)
INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR}/packages/epetra/src)

# Generate SWIG module
MAKE_SWIG(MODULE Epetra
    SOURCE Epetra.i
    DEPLIBS ${Trilinos_LIBRARIES}
)
