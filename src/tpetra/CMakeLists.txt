## Copyright 2017-2018, UT-Battelle, LLC
##
## SPDX-License-Identifier: BSD-3-Clause
## License-Filename: LICENSE

##---------------------------------------------------------------------------##
## Define the package
##---------------------------------------------------------------------------##

TRIBITS_SUBPACKAGE(Tpetra)

##---------------------------------------------------------------------------##
## Assert that Tpetra has ETI on the required types
##---------------------------------------------------------------------------##

# Make sure that Tpetra actually defined these variables, even if they are
# empty. This works even if ETI is off, as Tpetra defines these for running
# tests. When ETI is off, the ForTrilinos' required combination is on.
ASSERT_DEFINED(TpetraCore_ETI_SCALARS)
ASSERT_DEFINED(TpetraCore_ETI_LORDS)
ASSERT_DEFINED(TpetraCore_ETI_GORDS)
ASSERT_DEFINED(TpetraCore_ETI_NODES)

LIST (FIND TpetraCore_ETI_SCALARS "double" index)
IF (-1 EQUAL index) # not in the list of integer types
  MESSAGE(FATAL_ERROR "ForTrilinos requires Tpetra ETI for <double,int,long long,Kokkos::Compat::${Kokkos_NODE_TYPE}>")
ENDIF ()
LIST (FIND TpetraCore_ETI_LORDS "int" index)
IF (-1 EQUAL index) # not in the list of integer types
  MESSAGE(FATAL_ERROR "ForTrilinos requires Tpetra ETI for <double,int,long long,Kokkos::Compat::${Kokkos_NODE_TYPE}>")
ENDIF ()
LIST (FIND TpetraCore_ETI_GORDS "long long" index)
IF (-1 EQUAL index) # not in the list of integer types
  MESSAGE(FATAL_ERROR "ForTrilinos requires Tpetra ETI for <double,int,long long,Kokkos::Compat::${Kokkos_NODE_TYPE}>")
ENDIF ()

LIST (FIND TpetraCore_ETI_NODES "Kokkos::Compat::${Kokkos_NODE_TYPE}" index)
IF (-1 EQUAL index) # not in the list of integer types
  MESSAGE(FATAL_ERROR "ForTrilinos requires Tpetra ETI for <double,int,long long,Kokkos::Compat::${Kokkos_NODE_TYPE}>")
ENDIF ()

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

TRIBITS_ADD_TEST_DIRECTORIES(tutorials)

#TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

##---------------------------------------------------------------------------##
## D) Do standard postprocessing
##---------------------------------------------------------------------------##

TRIBITS_SUBPACKAGE_POSTPROCESS()
