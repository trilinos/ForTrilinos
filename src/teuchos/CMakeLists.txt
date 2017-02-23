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


IF (${PACKAGE_NAME}_ENABLE_DeveloperMode)
  # Generate wrapper and proxy files
  SET(CMAKE_SWIG_OUTDIR "${CMAKE_CURRENT_SOURCE_DIR}")

  MAKE_SWIG(MODULE forteuchos
    SOURCE forteuchos.i
    )

ELSE()
  # Use pre-generated files
  SET(SOURCES
    forteuchos.f90
    forteuchosFORTRAN_wrap.cxx)

  TRIBITS_ADD_LIBRARY(
    ForTeuchos.so
    SOURCES ${SOURCES})
ENDIF()

##---------------------------------------------------------------------------##
## end of teuchos/CMakeLists.txt
##---------------------------------------------------------------------------##
