##---------------------------------------------------------------------------##
## Define the package
##---------------------------------------------------------------------------##

TRIBITS_PACKAGE(ForTrilinos)

##---------------------------------------------------------------------------##
## Set up package-specific options
##---------------------------------------------------------------------------##

TRIBITS_ADD_DEBUG_OPTION()
TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

# Enable unlimited-length lines and C-like preprocessing
SET(CMAKE_Fortran_FLAGS "-ffree-line-length-none -cpp ${CMAKE_Fortran_FLAGS}")

IF (${PACKAGE_NAME}_ENABLE_DeveloperMode)
    # SWIG setup
    FIND_PACKAGE(SWIG REQUIRED)
    INCLUDE(SwigModules)
    # Ignore some SWIG warnings:
    #   401: "Nothing known about base class"
    SET(CMAKE_SWIG_FLAGS "-w401 ${CMAKE_SWIG_FLAGS}")

    IF (NOT TPL_ENABLE_MPI)
        # Warn if MPI is disabled since not all wrapper functions will be
        # generated. Later this might be extended to other optional packages.
        MESSAGE(WARNING "DeveloperMode is being used without MPI: "
          "the newly generated wrappers will NOT be generated in the source "
          "directory, since they will not contain MPI-related classes.")
    ELSE()
        # Export the newly generated wrappers to the source directory.
        SET(ForTrilinos_EXPORT_SWIG TRUE)
    ENDIF()
ENDIF()

##---------------------------------------------------------------------------##
## Add the libraries, tests, and examples
##---------------------------------------------------------------------------##

# TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)
# TRIBITS_ADD_TEST_DIRECTORIES(test)

##---------------------------------------------------------------------------##
## Build the documentation published on http://datatransferkit.readthedocs.org
##---------------------------------------------------------------------------##

IF(DataTransferKit_ENABLE_ReadTheDocs)
    ADD_SUBDIRECTORY(docs)
ENDIF()

##---------------------------------------------------------------------------##
## Generate the Doxygen documentation
##---------------------------------------------------------------------------##

IF(DataTransferKit_ENABLE_Doxygen)
    ADD_SUBDIRECTORY(docs/doxygen)
ENDIF()

TRIBITS_PROCESS_SUBPACKAGES()
TRIBITS_PACKAGE_DEF()
TRIBITS_PACKAGE_POSTPROCESS()

##---------------------------------------------------------------------------##
## Build the documentation published on http://fortrilinos.readthedocs.org
##---------------------------------------------------------------------------##
IF(ForTrilinos_ENABLE_ReadTheDocs)
    ADD_SUBDIRECTORY(docs)
ENDIF()

SET(ForTrilinos_PACKAGES "")
SET(ForTrilinos_MODULES  "")


#IF(ForTrilinos_ENABLE_Tpetra)
  APPEND_SET(ForTrilinos_PACKAGES Tpetra)
  APPEND_SET(ForTrilinos_MODULES  Tpetra)
#ENDIF(ForTrilinos_ENABLE_Tpetra)

#IF(ForTrilinos_ENABLE_Teuchos)
  APPEND_SET(ForTrilinos_PACKAGES Teuchos)
  APPEND_SET(ForTrilinos_MODULES  Teuchos)
#ENDIF(ForTrilinos_ENABLE_Teuchos)

# The default status for whether we build docstrings for ForTrilinos is
# based on whether the user has doxygen installed.
FIND_PACKAGE(Doxygen)
IF(DOXYGEN_FOUND)
  IF(${DOXYGEN_VERSION} VERSION_EQUAL 1.8.13 AND ForTrilinos_DOCSTRINGS)
    MESSAGE(FATAL_ERROR 
      "Doxygen version ${DOXYGEN_VERSION} has a bug that makes it incompatible with "
      "ForTrilinos. To work around this, you may disable ForTrilinos docstrings by setting "
      "ForTrilinos_DOCSTRINGS=OFF. "
      "Alternatively, you may choose to use a different version of doxygen. If you "
      "choose this route, be sure that the new version of doxygen is found first in "
      "your environment path.")
  ELSE(${DOXYGEN_VERSION} VERSION_EQUAL 1.8.13 AND ForTrilinos_DOCSTRINGS)
    SET(ForTrilinos_DOCSTRINGS_DEFAULT ON)
  ENDIF(${DOXYGEN_VERSION} VERSION_EQUAL 1.8.13 AND ForTrilinos_DOCSTRINGS)
ELSE(DOXYGEN_FOUND)
  SET(ForTrilinos_DOCSTRINGS_DEFAULT OFF)
  MESSAGE(STATUS "doxygen not found!")
ENDIF(DOXYGEN_FOUND)

# Set the ForTrilinos docstrings flag
SET(ForTrilinos_DOCSTRINGS ${ForTrilinos_DOCSTRINGS_DEFAULT} CACHE BOOL
    "a flag that determines whether Python docstrings will be built using doxygen")
IF(ForTrilinos_DOCSTRINGS AND NOT DOXYGEN_FOUND)
  MESSAGE(ERROR "ForTrilinos docstrings enabled but doxygen not found")
ENDIF(ForTrilinos_DOCSTRINGS AND NOT DOXYGEN_FOUND)

ADD_SUBDIRECTORY(doc)

