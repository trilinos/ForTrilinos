##---------------------------------------------------------------------------##
## File  : teuchos/CMakeLists.txt
## Author: Seth R Johnson
## Date  : Tue Feb 21 14:30:37 2017
## Brief : TriBITS package file for teuchos/
## Copyright (c) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
##---------------------------------------------------------------------------##

IF (${PACKAGE_NAME}_ENABLE_DeveloperMode)
  # Generate wrapper and proxy files
  SET(CMAKE_SWIG_OUTDIR "${CMAKE_CURRENT_SOURCE_DIR}")

  MAKE_SWIG(MODULE Teuchos
    SOURCE Teuchos.i
    DEPLIBS ${Trilinos_LIBRARIES}
    )
ELSE()
  # Use pre-generated files
  SET(SOURCES
    teuchos.f90
    teuchosFORTRAN_wrap.cxx)

  TRIBITS_ADD_LIBRARY(
    Teuchos.so
    SOURCES ${SOURCES}
    DEPLIBS ${Trilinos_LIBRARIES}
    )
ENDIF()

##---------------------------------------------------------------------------##
## end of teuchos/CMakeLists.txt
##---------------------------------------------------------------------------##
