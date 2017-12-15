!Copyright 2017, UT-Battelle, LLC
!
!SPDX-License-Identifier: BSD-3-Clause
!License-Filename: LICENSE
program test_TeuchosComm
#include "FortranTestUtilities.h"
#include "ForTrilinosTeuchos_config.hpp"
  use iso_fortran_env
  use, intrinsic :: iso_c_binding
  use forteuchos
#ifdef HAVE_MPI
    use mpi
#endif

  implicit none
  character(len=26), parameter :: FILENAME='test_teuchos_comm.f90'

  SETUP_TEST()

  ADD_SUBTEST_AND_RUN(TeuchosComm_Basic)

  TEARDOWN_TEST()
contains

  FORTRILINOS_UNIT_TEST(TeuchosComm_Basic)
    type(TeuchosComm) :: comm
    integer :: comm_rank_f, comm_rank_c
    integer :: comm_size_f, comm_size_c

#ifdef HAVE_MPI
    call comm%create(MPI_COMM_WORLD); TEST_IERR()
    TEST_EQUALITY_CONST(c_associated(comm%swigptr), .true.)

    call MPI_COMM_RANK(MPI_COMM_WORLD, comm_rank_f, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, comm_size_f, ierr)
#else
    call comm%create(); TEST_IERR()
    TEST_EQUALITY_CONST(c_associated(comm%swigptr), .true.)

    comm_rank_f = 0
    comm_size_f = 1
#endif

    comm_rank_c = comm%getRank(); TEST_IERR()
    TEST_EQUALITY(comm_rank_f, comm_rank_c)

    comm_size_c = comm%getSize(); TEST_IERR()
    TEST_EQUALITY(comm_size_f, comm_size_c)

    call comm%release()

  END_FORTRILINOS_UNIT_TEST(TeuchosComm_Basic)

end program test_TeuchosComm
