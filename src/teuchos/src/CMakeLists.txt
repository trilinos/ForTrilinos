##---------------------------------------------------------------------------##
## File  : teuchos/CMakeLists.txt
## Author: Seth R Johnson
## Date  : Tue Feb 21 14:30:37 2017
## Brief : TriBITS package file for teuchos/
## Copyright (c) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
##---------------------------------------------------------------------------##

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

IF (ForTrilinos_ENABLE_DeveloperMode)
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
    forteuchos
    SOURCES ${GEN_SOURCE} ${EXTRA_SOURCE})
ENDIF()

INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/forteuchos.mod DESTINATION include)

##---------------------------------------------------------------------------##
## end of teuchos/CMakeLists.txt
##---------------------------------------------------------------------------##
