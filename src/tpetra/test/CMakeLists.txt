## Copyright 2017-2018, UT-Battelle, LLC
##
## SPDX-License-Identifier: BSD-3-Clause
## License-Filename: LICENSE

##---------------------------------------------------------------------------##
## TESTS
##---------------------------------------------------------------------------##

SET(TESTS
  test_tpetra_map.F90
  test_tpetra_crsgraph.F90
  test_tpetra_crsmatrix.F90
  test_tpetra_multivector.F90
  test_tpetra_import_export.F90
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
