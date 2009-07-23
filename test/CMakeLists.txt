
INCLUDE(PackageAddExecutableAndTest)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  simpleEpetraOOFortranExample
  SOURCES simpleEpetraOOFortranExample.f90 interoperability_check.f90
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS -v
  STANDARD_PASS_OUTPUT
  )

#SET_SOURCE_FILES_PROPERTIES( ${SOURCES} PROPERTIES LANGUAGE Fortran )

SET_TARGET_PROPERTIES(ForTrilinos_simpleEpetraOOFortranExample PROPERTIES LINKER_LANGUAGE Fortran)

#PACKAGE_ADD_EXECUTABLE_AND_TEST(
#  interoperability_check
#  SOURCES interoperability_check.f90
#  COMM serial mpi
#  NUM_MPI_PROCS 1
#  ARGS -v
#  STANDARD_PASS_OUTPUT
#  )

#SET_SOURCE_FILES_PROPERTIES( ${SOURCES} PROPERTIES LANGUAGE Fortran )

#SET_TARGET_PROPERTIES(ForTrilinos_interoperability_check PROPERTIES LINKER_LANGUAGE Fortran)
