INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_INCLUDE_DIRS})

IF (ForTrilinos_ENABLE_DeveloperMode)
    # Generate wrapper and proxy files
    SET(CMAKE_SWIG_OUTDIR "${CMAKE_CURRENT_SOURCE_DIR}")

    MAKE_SWIG(MODULE forepetra
        SOURCE forepetra.i
        DEPLIBS ${Trilinos_LIBRARIES}
    )
ELSE()
    APPEND_GLOB(SOURCES *.f90)
    APPEND_GLOB(SOURCES *.cxx)

    TRIBITS_ADD_LIBRARY(
      forepetra
      SOURCES ${SOURCES}
      )
ENDIF()
