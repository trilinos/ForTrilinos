TRIBITS_COPY_FILES_TO_BINARY_DIR(Simple_cp
  SOURCE_FILES stratimikos.xml
 )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  trilinos_handle_test_cpp
  SOURCES trilinos_handle_test.cpp
  COMM serial mpi
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  trilinos_handle_test_fortran
  SOURCES trilinos_handle_test.f90 DBCF.f90
  COMM serial mpi
)

IF (${PACKAGE_NAME}_ENABLE_Belos AND ${PACKAGE_NAME}_ENABLE_Ifpack2)
  TRIBITS_ADD_TEST(
    trilinos_handle_test_cpp
    COMM serial mpi
    NUM_MPI_PROCS 1
  )

  TRIBITS_ADD_TEST(
    trilinos_handle_test_fortran
    COMM serial mpi
    NUM_MPI_PROCS 1
  )
ENDIF()
