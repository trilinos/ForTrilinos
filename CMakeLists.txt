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
