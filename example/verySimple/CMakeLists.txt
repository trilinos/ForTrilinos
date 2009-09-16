
INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  verySimple
  SOURCES verySimple.F90
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS -v
  STANDARD_PASS_OUTPUT
  )

#SET_SOURCE_FILES_PROPERTIES( ${SOURCES} PROPERTIES LANGUAGE Fortran )

# linker flag below only needed for gfortran!
SET_TARGET_PROPERTIES(ForTrilinos_verySimple PROPERTIES LINK_FLAGS -lgfortranbegin)
