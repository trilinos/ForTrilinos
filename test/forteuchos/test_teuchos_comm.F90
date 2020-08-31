!Copyright 2017-2018, UT-Battelle, LLC
!
!SPDX-License-Identifier: BSD-3-Clause
!License-Filename: LICENSE
program test_TeuchosComm
#include "FortranTestUtilities.h"
#include "ForTrilinos_config.h"
  use iso_fortran_env
  use, intrinsic :: iso_c_binding
#if FORTRILINOS_USE_MPI
    use mpi
#endif

#include "ForTrilinos.h"
  use forteuchos

  implicit none
  character(len=26), parameter :: FILENAME='test_teuchos_comm.F90'
  integer :: ierr, mpicomm


  SETUP_TEST()

  ADD_SUBTEST_AND_RUN(TeuchosComm_Basic)

  ! Unit tests, assume user has MPI enabled
  ADD_SUBTEST_AND_RUN(TeuchosComm_getRank)
  ADD_SUBTEST_AND_RUN(TeuchosComm_getSize)
  ADD_SUBTEST_AND_RUN(TeuchosComm_barrier)
  ADD_SUBTEST_AND_RUN(TeuchosComm_getRawMpiComm)

  TEARDOWN_TEST()
contains

  FORTRILINOS_UNIT_TEST(TeuchosComm_Basic)
    type(TeuchosComm) :: comm
    integer :: comm_rank_f, comm_rank_c
    integer :: comm_size_f, comm_size_c

#if FORTRILINOS_USE_MPI
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
#if FORTRILINOS_USE_MPI
    TEST_EQUALITY(MPI_COMM_WORLD, mpicomm)
#else
    ! Should have thrown an error
    TEST_INEQUALITY(0, fortrilinos_ierr)
    FORTRILINOS_IERR = 0
#endif

    call comm%release()

  END_FORTRILINOS_UNIT_TEST(TeuchosComm_Basic)

  ! ---------------------------------getRank---------------------------------- !
  FORTRILINOS_UNIT_TEST(TeuchosComm_getRank)
    type(TeuchosComm) :: comm
    integer :: myrank_t, myrank_m, ierr
    OUT0("Starting TeuchosComm_getRank!")

#if FORTRILINOS_USE_MPI
    comm = TeuchosComm(MPI_COMM_WORLD); TEST_IERR()
    myrank_t = comm%getRank(); TEST_IERR()
    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank_m, ierr)

    TEST_EQUALITY(myrank_t, myrank_m)

    call comm%release(); TEST_IERR()
#else
    return
#endif

    OUT0("Finished TeuchosComm_getRank!")

  END_FORTRILINOS_UNIT_TEST(TeuchosComm_getRank)

  ! ---------------------------------getSize---------------------------------- !
  FORTRILINOS_UNIT_TEST(TeuchosComm_getSize)
    type(TeuchosComm) :: comm
    integer :: mysize_t, mysize_m, ierr
    OUT0("Starting TeuchosComm_getSize!")

#if FORTRILINOS_USE_MPI
    comm = TeuchosComm(MPI_COMM_WORLD); TEST_IERR()
    mysize_t = comm%getSize(); TEST_IERR()
    call MPI_COMM_SIZE(MPI_COMM_WORLD, mysize_m, ierr)

    TEST_EQUALITY(mysize_t, mysize_m)

    call comm%release(); TEST_IERR()
#else
    return
#endif
    OUT0("Finished TeuchosComm_getSize!")

  END_FORTRILINOS_UNIT_TEST(TeuchosComm_getSize)

  ! ---------------------------------barrier---------------------------------- !
  FORTRILINOS_UNIT_TEST(TeuchosComm_barrier)
    type(TeuchosComm) :: comm
    OUT0("Starting TeuchosComm_barrier!")

#if FORTRILINOS_USE_MPI
    comm = TeuchosComm(MPI_COMM_WORLD); TEST_IERR()
    TEST_NOTHROW(call comm%barrier())

    call comm%release(); TEST_IERR()
#else
    return
#endif
    OUT0("Finished TeuchosComm_barrier!")

  END_FORTRILINOS_UNIT_TEST(TeuchosComm_barrier)

  ! ------------------------------getRawMpiComm------------------------------- !
  FORTRILINOS_UNIT_TEST(TeuchosComm_getRawMpiComm)
    type(TeuchosComm) :: comm
    OUT0("Starting TeuchosComm_getRawMpiComm!")

#if FORTRILINOS_USE_MPI
    comm = TeuchosComm(MPI_COMM_WORLD); TEST_IERR()
    TEST_EQUALITY(MPI_COMM_WORLD, comm%getRawMpiComm())

    call comm%release(); TEST_IERR()
#else
    return
#endif
    OUT0("Finished TeuchosComm_getRawMpiComm!")

  END_FORTRILINOS_UNIT_TEST(TeuchosComm_getRawMpiComm)

end program test_TeuchosComm
