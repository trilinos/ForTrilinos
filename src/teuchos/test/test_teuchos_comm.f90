!Copyright 2017-2018, UT-Battelle, LLC
!
!SPDX-License-Identifier: BSD-3-Clause
!License-Filename: LICENSE
program test_TeuchosComm
#include "FortranTestUtilities.h"
#include "ForTrilinosTeuchos_config.hpp"
  use iso_fortran_env
  use, intrinsic :: iso_c_binding
#ifdef HAVE_MPI
    use mpi
#endif

#include "ForTrilinos.h"
  use forteuchos

  implicit none
  character(len=26), parameter :: FILENAME='test_teuchos_comm.f90'
  integer :: ierr, mpicomm


  SETUP_TEST()

  ADD_SUBTEST_AND_RUN(TeuchosComm_Basic)

  TEARDOWN_TEST()
contains

  FORTRILINOS_UNIT_TEST(TeuchosComm_Basic)
    type(TeuchosComm) :: comm
    integer :: comm_rank_f, comm_rank_c
    integer :: comm_size_f, comm_size_c

#ifdef HAVE_MPI
    comm = TeuchosComm(MPI_COMM_WORLD); TEST_IERR()
    TEST_ASSERT(c_associated(comm%swigdata%cptr))

    call MPI_COMM_RANK(MPI_COMM_WORLD, comm_rank_f, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, comm_size_f, ierr)
#else
    comm = TeuchosComm(); TEST_IERR()
    TEST_ASSERT(c_associated(comm%swigdata%cptr))

    comm_rank_f = 0
    comm_size_f = 1
#endif

    comm_rank_c = comm%getRank(); TEST_IERR()
    TEST_EQUALITY(comm_rank_f, comm_rank_c)

    comm_size_c = comm%getSize(); TEST_IERR()
    TEST_EQUALITY(comm_size_f, comm_size_c)

    mpicomm = comm%getRawMpiComm()
#ifdef HAVE_MPI
    TEST_EQUALITY(MPI_COMM_WORLD, mpicomm)
#else
    ! Should have thrown an error
    TEST_INEQUALITY(0, fortrilinos_ierr)
    FORTRILINOS_IERR = 0
#endif

    call comm%release()

  END_FORTRILINOS_UNIT_TEST(TeuchosComm_Basic)

end program test_TeuchosComm
