##---------------------------------------------------------------------------##
## TESTS
##---------------------------------------------------------------------------##
SET(TESTS
  test_map.f90
  )

FOREACH(test ${TESTS})
    GET_FILENAME_COMPONENT(name ${test} NAME_WE)
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
        ${name}
        SOURCES ${test} DBCF.f90
        LINKER_LANGUAGE Fortran
        COMM serial mpi
        NUM_MPI_PROCS 4
    )
ENDFOREACH()
