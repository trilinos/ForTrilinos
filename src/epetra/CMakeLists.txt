INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_INCLUDE_DIRS})

IF (${PACKAGE_NAME}_ENABLE_DeveloperMode)
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
      ForEpetra.so
      SOURCES ${SOURCES}
      )
ENDIF()

