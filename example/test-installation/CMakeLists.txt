#----------------------------------*-CMake-*----------------------------------#
# Copyright 2020 UT-Battelle, LLC
# License-Filename: LICENSE
# SPDX-License-Identifier: BSD-3-Clause
#-----------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.12)
project(ForTrilinosInstallTest VERSION 0.0.1 LANGUAGES Fortran)

#-----------------------------------------------------------------------------#

find_package(ForTrilinos)

#-----------------------------------------------------------------------------#

add_executable(downstream-app downstream-app.F90)
target_link_libraries(downstream-app ForTrilinos::ForTrilinos MPI::MPI_Fortran)

include(CTest)
if(BUILD_TESTING)
  set(_EXE_NAME "$<TARGET_FILE:downstream-app>")
  if(ForTrilinos_USE_MPI)
    add_test(NAME "downstream-app" COMMAND
      "${MPIEXEC_EXECUTABLE}" ${MPIEXEC_NUMPROC_FLAG}
      "${ForTrilinos_MAX_NUMPROCS}"
       ${MPIEXEC_PREFLAGS} "${_EXE_NAME}" ${MPIEXEC_POSTFLAGS}
    )
  else()
    add_test(NAME "downstream-app" COMMAND "${_EXE_NAME}")
  endif()
endif()

#-----------------------------------------------------------------------------#

