## Copyright 2017, UT-Battelle, LLC
##
## SPDX-License-Identifier: BSD-3-Clause
## License-Filename: LICENSE
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

TRIBITS_ADD_LIBRARY(
  forutils
  HEADERS DBCF.h FortranTestMacros.h FortranTestUtilities.h
  SOURCES DBCF.f90 fortest.f90
)
