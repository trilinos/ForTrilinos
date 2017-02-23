INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/..)  # For ForTrilinosSimpleInterface_config.hpp

APPEND_GLOB(SOURCES *.cpp)
APPEND_GLOB(SOURCES *.cxx)
APPEND_GLOB(SOURCES *.f90)
APPEND_GLOB(HEADERS *.hpp)

TRIBITS_ADD_LIBRARY(
  simple_interface
  SOURCES simple_interface.cpp simple_interface_wrap.cxx simple_interface.f90
  )
