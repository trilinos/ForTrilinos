INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/../src)

TRIBITS_COPY_FILES_TO_BINARY_DIR(Simple_cp
  SOURCE_FILES stratimikos.xml
 )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  simple_interface_test_cpp
  SOURCES simple_interface_test.cpp
  COMM serial mpi
  DEPLIBS simple_interface
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  simple_interface_test_fortran
  SOURCES simple_interface_test.f90
  COMM serial mpi
  DEPLIBS simple_interface
  STANDARD_PASS_OUTPUT
)
