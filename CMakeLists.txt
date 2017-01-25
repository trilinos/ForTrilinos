
#
# A) Define the package
#

TRIBITS_PACKAGE(ForTrilinos CXX Fortran)

#
# B) Set up package-specific options
#
TRIBITS_ADD_DEBUG_OPTION()
TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()


SET(CMAKE_Fortran_FLAGS "-ffree-line-length-none ${CMAKE_Fortran_FLAGS}")

#
# C) Add the libraries, tests, and examples
#

# SWIG setup
FIND_PACKAGE(SWIG REQUIRED)
INCLUDE(SwigModules)

ADD_SUBDIRECTORY(src)

# TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)
# TRIBITS_ADD_TEST_DIRECTORIES(test)

#
# Exclude files for source package.
#

TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
#
#
