##---------------------------------------------------------------------------##
## File  : teuchos/CMakeLists.txt
## Author: Seth R Johnson
## Date  : Tue Feb 21 14:30:37 2017
## Brief : TriBITS package file for teuchos/
## Copyright (c) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
##---------------------------------------------------------------------------##

# Include the build directory so downstream Fortran code (including the unit
# tests) can use the generated modules
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
# Include the source dir to allow access to headers
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(EXTRA_SOURCE
  DBCF.f90)

IF (${PACKAGE_NAME}_ENABLE_DeveloperMode)
  # Generate wrapper and proxy files
  SET(CMAKE_SWIG_OUTDIR "${CMAKE_CURRENT_SOURCE_DIR}")

  MAKE_SWIG(MODULE forteuchos
    SOURCE forteuchos.i
    EXTRASRC ${EXTRA_SOURCE}
    )

ELSE()
  # Use pre-generated files
  SET(GEN_SOURCE
    forteuchos.f90
    forteuchosFORTRAN_wrap.cxx)

  TRIBITS_ADD_LIBRARY(
    ForTeuchos.so
    SOURCES ${GEN_SOURCE} ${EXTRA_SOURCE})
ENDIF()

TRIBITS_ADD_TEST_DIRECTORIES(test)

##---------------------------------------------------------------------------##
## end of teuchos/CMakeLists.txt
##---------------------------------------------------------------------------##
