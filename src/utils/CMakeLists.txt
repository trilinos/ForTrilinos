TRIBITS_SUBPACKAGE(Utils)

TRIBITS_CONFIGURE_FILE(FortranTestMacros.h)

##---------------------------------------------------------------------------##
## Generate config file and add it to the list
##---------------------------------------------------------------------------##
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

##---------------------------------------------------------------------------##
## Add library, test, and examples.
##---------------------------------------------------------------------------##

ADD_SUBDIRECTORY(src)

# TRIBITS_ADD_TEST_DIRECTORIES(test)

#TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

##---------------------------------------------------------------------------##
## D) Do standard postprocessing
##---------------------------------------------------------------------------##

TRIBITS_SUBPACKAGE_POSTPROCESS()
