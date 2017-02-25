TRIBITS_COPY_FILES_TO_BINARY_DIR(Simple_cp
  SOURCE_FILES stratimikos.xml
 )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  simple_interface_test_cpp
  SOURCES simple_interface_test.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  simple_interface_test_fortran
  SOURCES simple_interface_test.f90 DBCF.f90
  COMM serial mpi
  NUM_MPI_PROCS 1
)
