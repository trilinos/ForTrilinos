##---------------------------------------------------------------------------##
## File  : teuchos/CMakeLists.txt
## Author: Seth R Johnson
## Date  : Tue Feb 21 14:30:37 2017
## Brief : TriBITS package file for teuchos/
## Copyright (c) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
##---------------------------------------------------------------------------##

TRIBITS_SUBPACKAGE(Belos)

##---------------------------------------------------------------------------##
## Generate config file and add it to the list
##---------------------------------------------------------------------------##
TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.hpp)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

##---------------------------------------------------------------------------##
## Add library, test, and examples.
##---------------------------------------------------------------------------##

ADD_SUBDIRECTORY(src)

TRIBITS_ADD_TEST_DIRECTORIES(test)

#TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

##---------------------------------------------------------------------------##
## D) Do standard postprocessing
##---------------------------------------------------------------------------##

TRIBITS_SUBPACKAGE_POSTPROCESS()
