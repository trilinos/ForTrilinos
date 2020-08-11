#---------------------------------*-CMake-*----------------------------------#
# Copyright 2020 UT-Battelle, LLC and ForTrilinos developers.
# SPDX-License-Identifier: BSD-3-Clause
# License-Filename: LICENSE
#[=======================================================================[.rst:

FindTrilinos
------------

Find Trilinos and define modern CMake targets.

.. code-block:: cmake

  find_package(Trilinos REQUIRED COMPONENTS TeuchosNumerics)
  target_link_libraries(<MYTARGET> Trilinos::TeuchosNumerics)

This script changes Trilinos targets from library-name-based to component-based.

#]=======================================================================]

include(FindPackageHandleStandardArgs)

# Load any installed $TRILINOS_ROOT/lib/cmake/Trilinos/TrilinosConfig.cmake
find_package(Trilinos QUIET NO_MODULE)

if(Trilinos_FOUND AND NOT TARGET "Trilinos::Trilinos")
  # Add base-level trilinos includes
  add_library(Trilinos::Trilinos INTERFACE IMPORTED)
  target_include_directories(Trilinos::Trilinos SYSTEM INTERFACE
    "${Trilinos_INCLUDE_DIRS}"
  )

  set(_components "${Trilinos_FIND_COMPONENTS}")
  if(NOT _components)
    # Set to variable leaked by TrilinosConfig.cmake with all top-level
    # packages.
    set(_components "${COMPONENTS_LIST}")
  endif()

  foreach(_comp IN LISTS _components)
    if(NOT Trilinos_${_comp}_FOUND)
      continue()
    endif()
    # We want to provide targets based on component names, but Trilinos only
    # makes targets with the lowercase library names. We can't use ALIAS because
    # Trilinos also fails to declare the targets as GLOBAL. A solution is to
    # declare our own scoped library and use interface dependencies to link
    # against the first library name for that component.
    add_library("Trilinos::${_comp}" INTERFACE IMPORTED)

    # Link library against the Trilinos-created target and the main Trilinos
    # target (for include dirs/flags)
    list(GET ${_comp}_LIBRARIES 0 _lib)
    target_link_libraries("Trilinos::${_comp}"
      INTERFACE Trilinos::Trilinos "${_lib}")
  endforeach()
  unset(_components)
  unset(_comp)
  unset(_lib)
endif()


if(Trilinos_Tpetra_FOUND AND PROJECT_NAME STREQUAL "ForTrilinos")
  set(_require_tpetra_inst Tpetra_INST_INT_LONG_LONG)
  set(_tpetra_check "${Trilinos_DIR}")
  if(NOT Tpetra_INST_INT_LONG_LONG STREQUAL "${_tpetra_check}")
    try_compile(Tpetra_INST_INT_LONG_LONG
      "${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}"
      "${PROJECT_SOURCE_DIR}/cmake/FindTrilinos/CheckTpetraInstantiations.cc"
      LINK_LIBRARIES Trilinos::Tpetra
      OUTPUT_VARIABLE _tpetra_out
    )
    if(NOT Tpetra_INST_INT_LONG_LONG)
      set(_tpetra_check "FALSE")
      message(STATUS "ForTrilinos requires Trilinos built "
        "with <LO=int, GO=long long>: ${_tpetra_out}")
    endif()
    set(Tpetra_INST_INT_LONG_LONG "${_tpetra_check}" CACHE INTERNAL
      "Verified path to Trilinos")
    unset(_tpetra_check)
    unset(_tpetra_out)
  endif()
endif()

find_package_handle_standard_args(
  Trilinos HANDLE_COMPONENTS
  REQUIRED_VARS Trilinos_INCLUDE_DIRS Trilinos_LIBRARY_DIRS Trilinos_LIBRARIES
    ${_require_tpetra_inst}
  VERSION_VAR Trilinos_VERSION
)

unset(_require_tpetra_inst)

#----------------------------------------------------------------------------#
