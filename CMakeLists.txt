## Copyright 2017-2018, UT-Battelle, LLC
##
## SPDX-License-Identifier: BSD-3-Clause
## License-Filename: LICENSE

##---------------------------------------------------------------------------##
## Define the package
##---------------------------------------------------------------------------##

TRIBITS_PACKAGE_DECL(ForTrilinos)

##---------------------------------------------------------------------------##
## Set up package-specific options
##---------------------------------------------------------------------------##

TRIBITS_ADD_DEBUG_OPTION()
TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

# Pick a default node to compile for in order of Cuda -> OpenMP -> Serial
# For now we just build for this one node which is selected based on the config.
# This is an intermediate step and expected to change later for supporting
# all nodes simultaneously.
IF (KOKKOS_ENABLE_CUDA)
  SET(Kokkos_NODE_TYPE "KokkosCudaWrapperNode")
ELSEIF (KOKKOS_ENABLE_OPENMP)
  SET(Kokkos_NODE_TYPE "KokkosOpenMPWrapperNode")
ELSE ()
  SET(Kokkos_NODE_TYPE "KokkosSerialWrapperNode")
ENDIF()
MESSAGE(STATUS "ForTrilinos default node type: ${Kokkos_NODE_TYPE}")

# Enable unlimited-length lines
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  SET(CMAKE_Fortran_FLAGS "-ffree-line-length-none ${CMAKE_Fortran_FLAGS}")
endif()

IF (${PACKAGE_NAME}_ENABLE_DeveloperMode)
    # SWIG setup
    FIND_PACKAGE(SWIG REQUIRED)
    # SWIG is requested and available; make sure it's the Fortran fork.
    INCLUDE(CheckSWIGFortran)
    IF (CMAKE_VERSION VERSION_LESS 3.10)
      # TODO: Old version of cmake modules
      include(SwigModulesOld)
    ELSE()
      include(SwigModules)
    ENDIF()

    # Ignore some SWIG warnings:
    #   401: "Nothing known about base class"
    LIST(APPEND CMAKE_SWIG_FLAGS "-w401" "-fext" "F90")
    set(SWIG_FORTRAN_EXTRA_FILE_EXTENSIONS ".F90")

    IF (NOT TPL_ENABLE_MPI)
        # Warn if MPI is disabled since not all wrapper functions will be
        # generated. Later this might be extended to other optional packages.
        MESSAGE(WARNING "DeveloperMode is being used without MPI: "
          "the newly generated wrappers will NOT be generated in the source "
          "directory, since they will not contain MPI-related classes.")
    ELSE()
        # Export the newly generated wrappers to the source directory.
        SET(ForTrilinos_EXPORT_SWIG TRUE CACHE BOOL
          "Build SWIG wrapper files in the source directory")
    ENDIF()
ENDIF()

##---------------------------------------------------------------------------##
## Add the libraries, tests, and examples
##---------------------------------------------------------------------------##

# TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)
# TRIBITS_ADD_TEST_DIRECTORIES(test)

##---------------------------------------------------------------------------##
## Generate the Doxygen documentation
##---------------------------------------------------------------------------##
# There does not currently exist Doxygen documentation
#IF(ForTrilinos_ENABLE_Doxygen)
#    ADD_SUBDIRECTORY(docs/doxygen)
#ENDIF()

##---------------------------------------------------------------------------##
## Build the documentation published on http://fortrilinos.readthedocs.org
##---------------------------------------------------------------------------##
IF(ForTrilinos_ENABLE_ReadTheDocs)
    ADD_SUBDIRECTORY(docs)
ENDIF()

TRIBITS_PROCESS_SUBPACKAGES()
TRIBITS_PACKAGE_DEF()
TRIBITS_PACKAGE_POSTPROCESS()
