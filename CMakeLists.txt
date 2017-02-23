
#
# A) Define the package
#

TRIBITS_PACKAGE(ForTrilinos CXX Fortran)

#
# B) Set up package-specific options
#
TRIBITS_ADD_DEBUG_OPTION()
TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()


# Enable unlimited-length lines and C-like preprocessing
SET(CMAKE_Fortran_FLAGS "-ffree-line-length-none -cpp ${CMAKE_Fortran_FLAGS}")

#
# C) Add the libraries, tests, and examples
#
IF (${PACKAGE_NAME}_ENABLE_DeveloperMode)
    # SWIG setup
    FIND_PACKAGE(SWIG REQUIRED)
    INCLUDE(SwigModules)
ENDIF()

# TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)
# TRIBITS_ADD_TEST_DIRECTORIES(test)

#
# Exclude files for source package.
#

TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

#
# D) Do standard postprocessing
#

TRIBITS_PROCESS_SUBPACKAGES()
TRIBITS_PACKAGE_DEF()
TRIBITS_PACKAGE_POSTPROCESS()
#
#
