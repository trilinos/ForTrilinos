!Copyright 2017-2018, UT-Battelle, LLC
!
!SPDX-License-Identifier: BSD-3-Clause
!License-Filename: LICENSE
program test_TrilinosSolver
#include "FortranTestUtilities.h"
  use iso_fortran_env
  use, intrinsic :: iso_c_binding
  use forteuchos
  use fortrilinos_hl
  use fortpetra
  use test_tpetra_crsmatrix_helper

  implicit none
  type(TeuchosComm) :: comm
  character(len=25), parameter :: FILENAME="test_tri_solve.F90"
  integer, parameter :: dp = kind(0.d0)
  integer :: my_rank, num_procs

  SETUP_TEST()

#if FORTRILINOS_USE_MPI
  comm = TeuchosComm(MPI_COMM_WORLD); FORTRILINOS_CHECK_IERR()
#else
  comm = TeuchosComm()
#endif

  num_procs = comm%getSize()
  my_rank = comm%getRank()

  ADD_SUBTEST_AND_RUN(TrilinosSolver_setup_matrix)
  ADD_SUBTEST_AND_RUN(TrilinosSolver_setup_solver)
  ADD_SUBTEST_AND_RUN(TrilinosSolver_finalize)

  ! Abstract class - -  no need to test
  !ADD_SUBTEST_AND_RUN(TrilinosSolver_setup_operator)

  ADD_SUBTEST_AND_RUN(TrilinosSolver_solve)

  call comm%release()

  TEARDOWN_TEST()

contains

  ! -------------------------------setup_matrix------------------------------- !
  FORTRILINOS_UNIT_TEST(TrilinosSolver_setup_matrix)
    type(TrilinosSolver) :: eig
    type(TpetraMap) :: map
    type(TpetraCrsMatrix) :: a
    OUT0("Starting TrilinosSolver_setup_matrix!")

    call Tpetra_CrsMatrix_CreateIdentity(comm, a)
    eig = TrilinosSolver(); TEST_IERR()
    call eig%init(comm); TEST_IERR()
    TEST_NOTHROW(call eig%setup_matrix(a))

    call a%release(); TEST_IERR()
    call eig%finalize(); TEST_IERR()
    call eig%release(); TEST_IERR()

    OUT0("Finished TrilinosSolver_setup_matrix!")

  END_FORTRILINOS_UNIT_TEST(TrilinosSolver_setup_matrix)

  ! -------------------------------setup_solver------------------------------- !
  FORTRILINOS_UNIT_TEST(TrilinosSolver_setup_solver)
    type(TrilinosSolver) :: eig
    type(ParameterList) :: plist
    type(TpetraMap) :: map
    type(TpetraCrsMatrix) :: A
    OUT0("Starting TrilinosSolver_setup_solver!")

    call Tpetra_CrsMatrix_CreateIdentity(comm, A)
    call A%fillComplete(); FORTRILINOS_CHECK_IERR()

    plist = ParameterList("Stratimikos"); FORTRILINOS_CHECK_IERR()
    call load_from_xml(plist, 'stratimikos.xml')

    eig = TrilinosSolver(); TEST_IERR()
    call eig%init(comm); TEST_IERR()
    call eig%setup_matrix(A); TEST_IERR()

    TEST_NOTHROW(call eig%setup_solver(plist))

    call plist%release(); TEST_IERR()
    call A%release(); TEST_IERR()
    call eig%finalize(); TEST_IERR()
    call eig%release(); TEST_IERR()

    OUT0("Finished TrilinosSolver_setup_solver!")

  END_FORTRILINOS_UNIT_TEST(TrilinosSolver_setup_solver)
  ! ----------------------------------solve----------------------------------- !
  FORTRILINOS_UNIT_TEST(TrilinosSolver_solve)
    ! Yes, we do need to lift the beefy part of the previous test to test this function

    integer(global_size_type) :: n_global
    integer(size_type) :: num_vecs = 1

    integer :: ierr
    integer :: i, n
    integer(global_ordinal_type) :: offset

    type(ParameterList) :: plist
    type(TrilinosSolver) :: solver
    type(TpetraCrsMatrix) :: A
    type(TpetraMultiVector) :: B, X, residual
    type(TpetraMap) :: map

    real(scalar_type), dimension(:), allocatable :: rhs
    real(norm_type), dimension(:), allocatable :: norms
    real(scalar_type) :: r0, sone = 1., szero = 0.
    real(scalar_type), parameter :: tol=1.e-4

    OUT0("Starting TrilinosSolver_solve!")

    ! Read in the parameterList
    plist = ParameterList("Stratimikos"); FORTRILINOS_CHECK_IERR()
    call load_from_xml(plist, 'stratimikos.xml'); FORTRILINOS_CHECK_IERR()

    call  Tpetra_CrsMatrix_CreateTestMatrix_A(comm,A)
    call A%fillComplete(); FORTRILINOS_CHECK_IERR()
    n = A%getNodeNumRows()
    map = A%getRowMap()

    B = TpetraMultiVector(map, num_vecs); FORTRILINOS_CHECK_IERR()
    X = TpetraMultiVector(map, num_vecs); FORTRILINOS_CHECK_IERR()
    residual = TpetraMultiVector(map, num_vecs, .false.); FORTRILINOS_CHECK_IERR()

    allocate(norms(1))

    solver = TrilinosSolver(); FORTRILINOS_CHECK_IERR()

    call solver%init(comm); FORTRILINOS_CHECK_IERR()

    call solver%setup_matrix(A); FORTRILINOS_CHECK_IERR()

    call solver%setup_solver(plist); FORTRILINOS_CHECK_IERR()

    call X%randomize()
    call B%randomize()

    call A%apply(X, residual, TeuchosNO_TRANS, sone, szero); FORTRILINOS_CHECK_IERR()
    call residual%update(sone, B, -sone); FORTRILINOS_CHECK_IERR()
    call residual%norm2(norms); FORTRILINOS_CHECK_IERR()
    r0 = norms(1)

    ! Solve the system
    TEST_NOTHROW(call solver%solve(B, X))

     ! Check the solution
    call A%apply(X, residual, TeuchosNO_TRANS, sone, szero); FORTRILINOS_CHECK_IERR()
    call residual%update(sone, B, -sone); FORTRILINOS_CHECK_IERR()
    call residual%norm2(norms); FORTRILINOS_CHECK_IERR()

    if (norms(1)/r0 > tol) then
      write(error_unit, '(A)') 'The solver did not converge to the specified residual!'
      stop 1
    end if

   call solver%finalize(); FORTRILINOS_CHECK_IERR()

    call map%release(); FORTRILINOS_CHECK_IERR()
    call A%release(); FORTRILINOS_CHECK_IERR()
    call solver%release(); FORTRILINOS_CHECK_IERR()
    call plist%release(); FORTRILINOS_CHECK_IERR()
    call X%release(); FORTRILINOS_CHECK_IERR()
    call B%release(); FORTRILINOS_CHECK_IERR()
    call residual%release(); FORTRILINOS_CHECK_IERR()
    deallocate(norms)

    OUT0("Finished TrilinosSolver_solve!")

  END_FORTRILINOS_UNIT_TEST(TrilinosSolver_solve)
 ! ---------------------------------finalize--------------------------------- !
  FORTRILINOS_UNIT_TEST(TrilinosSolver_finalize)
    type(TrilinosSolver) :: eig
    OUT0("Starting TrilinosSolver_finalize!")

    eig = TrilinosSolver(); TEST_IERR()
    call eig%init(comm); TEST_IERR()
    TEST_NOTHROW(call eig%finalize())

    call eig%release(); TEST_IERR()

    OUT0("Finished TrilinosSolver_finalize!")

  END_FORTRILINOS_UNIT_TEST(TrilinosSolver_finalize)

end program test_TrilinosSolver
