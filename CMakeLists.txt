##---------------------------------------------------------------------------##
## Define the package
##---------------------------------------------------------------------------##

TRIBITS_PACKAGE(ForTrilinos CXX Fortran)

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
