INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/packages/epetra/src)
INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR}/packages/epetra/src)

# Generate SWIG module
MAKE_SWIG(MODULE Epetra0
    SOURCE Epetra_Object.i
    DEPLIBS ${Trilinos_LIBRARIES}
)
MAKE_SWIG(MODULE Epetra1
    SOURCE Epetra_Comm.i
    DEPLIBS ${Trilinos_LIBRARIES}
)
MAKE_SWIG(MODULE Epetra2
    SOURCE Epetra_SerialComm.i
    DEPLIBS ${Trilinos_LIBRARIES}
)
