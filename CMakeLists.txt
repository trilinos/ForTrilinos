#---------------------------------*-CMake-*----------------------------------#
# Copyright 2020 UT-Battelle, LLC
# License-Filename: LICENSE
# SPDX-License-Identifier: BSD-3-Clause
#----------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.12)
project(ForTrilinos VERSION 2.0.0 LANGUAGES Fortran CXX)
cmake_policy(VERSION 3.12...3.18)

list(APPEND CMAKE_MODULE_PATH
  "${CMAKE_CURRENT_LIST_DIR}/cmake"
)
if(CMAKE_VERSION VERSION_LESS 3.18)
  list(APPEND CMAKE_MODULE_PATH
    "${CMAKE_CURRENT_LIST_DIR}/cmake/backport/3.18"
  )
endif()

include(GNUInstallDirs)

#-----------------------------------------------------------------------------#
# Component options
#-----------------------------------------------------------------------------#

# Only enable SWIG fortran when Fortran is enabled
option(ForTrilinos_USE_SWIG_Fortran
  "Regenerate Fortran bindings with SWIG" OFF)

#-----------------------------------------------------------------------------#
# Deployment options
#-----------------------------------------------------------------------------#

option(BUILD_SHARED_LIBS "Default to building shared libraries" ON)

if(NOT DEFINED CMAKE_INSTALL_RPATH_USE_LINK_PATH)
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON CACHE STRING
    "Inform installed binaries of external library rpaths")
endif()
if(BUILD_SHARED_LIBS)
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}" CACHE STRING
    "Inform installed binaries of internal library rpaths")
endif()
if(APPLE)
  option(CMAKE_MACOSX_RPATH "Support @rpath in install targets" ON)
endif()

#-----------------------------------------------------------------------------#
# Test options
#-----------------------------------------------------------------------------#

option(ForTrilinos_TESTING "Build ForTrilinos unit tests" ON)
option(ForTrilinos_EXAMPLES "Build ForTrilinos examples" ON)
set(ForTrilinos_MAX_NUMPROCS "0" CACHE STRING
  "Maximum number of processors to use in tests")

#-----------------------------------------------------------------------------#
# Internal variables
#-----------------------------------------------------------------------------#

# Prefix for exported ForTrilinos targets
set(ForTrilinos_NAMESPACE "ForTrilinos::")

# Destination for configured header files (config.h) during build
set(ForTrilinos_HEADER_CONFIG_DIRECTORY "${PROJECT_BINARY_DIR}/include")
# Where to configure cmake files
set(ForTrilinos_CMAKE_CONFIG_DIRECTORY "${PROJECT_BINARY_DIR}/share")
# Where to build .mod files
set(ForTrilinos_Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/fortran")

# Where to install cmake config files
set(ForTrilinos_INSTALL_CMAKECONFIGDIR "${CMAKE_INSTALL_LIBDIR}/cmake/ForTrilinos")

#-----------------------------------------------------------------------------#
# LANGUAGE SUPPORT
#-----------------------------------------------------------------------------#

# Don't use compiler language extensions
set(CMAKE_CXX_EXTENSIONS OFF)

#-----------------------------------------------------------------------------#
# DEPENDENCIES
#
# Derive options from Trilinos configuration where possible.
#-----------------------------------------------------------------------------#

find_package(Trilinos 12.17 REQUIRED MODULE
  COMPONENTS
    TeuchosCore TeuchosComm TeuchosNumerics TeuchosParameterList
  OPTIONAL_COMPONENTS
    # Optional
    Belos
    Tpetra
    # Required for interface package
    Anasazi NOX Stratimikos Thyra ThyraTpetraAdapters
    # Optional for "interface" package
    Amesos2 Ifpack2 MueLu
)

set(Trilinos_COMPONENTS ${COMPONENTS_LIST})
if("MPI" IN_LIST Trilinos_TPL_LIST)
  set(Trilinos_USE_MPI ON)
endif()

set(ForTrilinos_USE_MPI "${Trilinos_USE_MPI}" CACHE BOOL
  "MPI status depends on the Trilinos installation in use" FORCE)

set(_fortrilinos_use_hl OFF)
if(Trilinos_Anasazi_FOUND
    AND Trilinos_NOX_FOUND
    AND Trilinos_Stratimikos_FOUND
    AND Trilinos_Thyra_FOUND
    AND Trilinos_ThyraTpetraAdapters_FOUND)
  set(_fortrilinos_use_hl ON)
endif()
set(ForTrilinos_USE_HL "${_fortrilinos_use_hl}" CACHE BOOL
  "High-level ForTrilinos wrappers depend on Trilinos installation" FORCE)

if(ForTrilinos_USE_MPI)
  if(NOT MPIEXEC_EXECUTABLE)
    # Hint for MPI
    set(MPIEXEC_EXECUTABLE "${TeuchosComm_MPI_EXEC}")
  endif()
  find_package(MPI REQUIRED COMPONENTS CXX Fortran)
  if(NOT ForTrilinos_MAX_NUMPROCS)
    # Not set or default of zero
    set(ForTrilinos_MAX_NUMPROCS "${MPIEXEC_MAX_NUMPROCS}" CACHE STRING
      "Maximum number of processors to use in tests" FORCE)
  elseif(ForTrilinos_MAX_NUMPROCS GREATER MPIEXEC_MAX_NUMPROCS)
    message(WARNING "Maximim number of processors in MPI "
      "tests exceeds the system-detected processor count: "
      "${ForTrilinos_MAX_NUMPROCS} > ${MPIEXEC_MAX_NUMPROCS}"
    )
  endif()
else()
  set(ForTrilinos_MAX_NUMPROCS "1" CACHE INTERNAL
    "Maximum number of processors to use in tests")
endif()

if(ForTrilinos_USE_SWIG_Fortran)
  find_package(SWIG REQUIRED COMPONENTS fortran)
  include(UseSWIG)

  if(NOT ForTrilinos_USE_MPI)
    message(WARNING "DeveloperMode is being used without MPI: "
      "the newly generated wrappers will NOT be generated in the source "
      "directory, since they will not contain MPI-related classes.")
  endif()
endif()

if(ForTrilinos_TESTING OR ForTrilinos_EXAMPLES)
  include(CTest)
endif()

#---------------------------------------------------------------------------##
# LIBRARY
#---------------------------------------------------------------------------##

add_subdirectory(src)

#---------------------------------------------------------------------------##
# TESTS
#---------------------------------------------------------------------------##

if(ForTrilinos_TESTING)
  add_subdirectory(test)
endif()

#---------------------------------------------------------------------------##
# EXAMPLES
#---------------------------------------------------------------------------##

if(ForTrilinos_EXAMPLES)
  add_subdirectory(example)
endif()

#---------------------------------------------------------------------------##
# INSTALLATION
#---------------------------------------------------------------------------##

include(ForTrilinosConfigSetup)
fortrilinos_configure_export()

#---------------------------------------------------------------------------##
