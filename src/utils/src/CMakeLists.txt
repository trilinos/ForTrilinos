## Copyright 2017-2018, UT-Battelle, LLC
##
## SPDX-License-Identifier: BSD-3-Clause
## License-Filename: LICENSE
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

IF (ForTrilinos_ENABLE_DeveloperMode)
  IF (ForTrilinos_EXPORT_SWIG)
    # Generate wrapper and proxy files since MPI is enabled
    SET(CMAKE_SWIG_OUTDIR "${CMAKE_CURRENT_SOURCE_DIR}")
  ENDIF()

  MAKE_SWIG(MODULE forerror
    SOURCE forerror.i
    EXTRASRC ${EXTRA_SOURCE}
    )
ELSE()
  # Use pre-generated files
  SET(GEN_SOURCE
    forerror.f90
    forerrorFORTRAN_wrap.cxx)

  TRIBITS_ADD_LIBRARY(
    forerror
    SOURCES ${GEN_SOURCE})
ENDIF()

TRIBITS_ADD_LIBRARY(
  forutils
  NOINSTALLHEADERS FortranTestUtilities.h
  SOURCES fortest.f90
  DEPLIBS forerror
)

INSTALL(FILES ${CMAKE_CURRENT_SOURCE_DIR}/ForTrilinos.h DESTINATION include)
INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/forerror.mod DESTINATION include)
