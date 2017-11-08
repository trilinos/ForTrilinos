##---------------------------------------------------------------------------##
## File  : teuchos/test/CMakeLists.txt
## Author: Seth R Johnson
## Date  : Thu Feb 23 10:47:45 2017
## Brief : TriBITS package file for teuchos/test
## Copyright (c) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
##---------------------------------------------------------------------------##

##---------------------------------------------------------------------------##
## TESTS
##---------------------------------------------------------------------------##

SET(TESTS
  test_array.f90
  test_arrayview.f90
  test_comm.f90
  test_plist.f90
  )

FOREACH(test ${TESTS})
    GET_FILENAME_COMPONENT(name ${test} NAME_WE)
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
        ${name}
        SOURCES ${test}
        LINKER_LANGUAGE Fortran
        COMM serial mpi
        NUM_MPI_PROCS 4
    )
ENDFOREACH()
