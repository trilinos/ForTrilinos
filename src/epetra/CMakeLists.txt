INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/packages/epetra/src)
INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR}/packages/epetra/src)

IF (${PACKAGE_NAME}_ENABLE_DeveloperMode)
    # Generate wrapper and proxy files
    SET(CMAKE_SWIG_OUTDIR "${CMAKE_CURRENT_SOURCE_DIR}")

    MAKE_SWIG(MODULE Epetra
        SOURCE Epetra.i
        DEPLIBS ${Trilinos_LIBRARIES}
    )
ELSE()
    APPEND_GLOB(SOURCES *.f90)
    APPEND_GLOB(SOURCES *.cxx)

    TRIBITS_ADD_LIBRARY(
      Epetra.so
      SOURCES ${SOURCES}
      DEPLIBS ${Trilinos_LIBRARIES}
      )
ENDIF()

