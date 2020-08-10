## Copyright 2017-2018, UT-Battelle, LLC
##
## SPDX-License-Identifier: BSD-3-Clause
## License-Filename: LICENSE

##---------------------------------------------------------------------------##
## TESTS
##---------------------------------------------------------------------------##

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(TESTS
  test_tpetra_map.F90
  test_tpetra_crsgraph.F90
  test_tpetra_crsmatrix.F90
  test_tpetra_multivector.F90
  test_tpetra_import_export.F90
  )
SET(TEST_HELPERS
  test_tpetra_crsgraph_helper.F90
  test_tpetra_crsmatrix_helper.F90
  test_tpetra_import_export_helper.F90
  test_tpetra_multivector_helper.F90
  )

TRIBITS_ADD_LIBRARY(
  fortpetra_helpers
  SOURCES ${TEST_HELPERS}
  DEPLIBS fortpetra
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
