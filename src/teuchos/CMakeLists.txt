##---------------------------------------------------------------------------##
## File  : teuchos/CMakeLists.txt
## Author: Seth R Johnson
## Date  : Tue Feb 21 14:30:37 2017
## Brief : TriBITS package file for teuchos/
## Copyright (c) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
##---------------------------------------------------------------------------##

##---------------------------------------------------------------------------##
## Define the package
##---------------------------------------------------------------------------##

TRIBITS_SUBPACKAGE(Teuchos)

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
