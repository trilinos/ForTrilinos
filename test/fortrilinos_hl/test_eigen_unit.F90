!Copyright 2017-2018, UT-Battelle, LLC
!
!SPDX-License-Identifier: BSD-3-Clause
!License-Filename: LICENSE
program test_TrilinosEigenSolver

#include "ForTrilinos_config.h"
#include "ForTrilinos.h"
#include "FortranTestUtilities.h"
  use iso_fortran_env
  use, intrinsic :: iso_c_binding
  use forteuchos
  use fortrilinos_hl
  use fortpetra
  use test_tpetra_crsmatrix_helper

  implicit none
  type(TeuchosComm) :: comm
  character(len=30), parameter :: FILENAME="test_eigen_unit.F90"
  integer, parameter :: dp = kind(0.d0)

  SETUP_TEST()

#if FORTRILINOS_USE_MPI
  comm = TeuchosComm(MPI_COMM_WORLD); FORTRILINOS_CHECK_IERR()
#else
  comm = TeuchosComm(); FORTRILINOS_CHECK_IERR()
#endif

  ADD_SUBTEST_AND_RUN(TrilinosEigenSolver_setup_matrix)
  ADD_SUBTEST_AND_RUN(TrilinosEigenSolver_setup_matrix_rhs)
  ADD_SUBTEST_AND_RUN(TrilinosEigenSolver_setup_solver)
  ADD_SUBTEST_AND_RUN(TrilinosEigenSolver_solve)

  ! No need to test at this level
  !ADD_SUBTEST_AND_RUN(TrilinosEigenSolver_setup_operator)
  !ADD_SUBTEST_AND_RUN(TrilinosEigenSolver_setup_operator_rhs)


  call comm%release()

  TEARDOWN_TEST()

contains

  ! -------------------------------setup_matrix------------------------------- !
  FORTRILINOS_UNIT_TEST(TrilinosEigenSolver_setup_matrix)
    type(TrilinosEigenSolver) :: eig
    type(TpetraCrsMatrix) :: a
    OUT0("Starting TrilinosEigenSolver_setup_matrix!")

    call Tpetra_CrsMatrix_CreateIdentity(comm, a)
    eig = TrilinosEigenSolver(comm); TEST_IERR()
    TEST_NOTHROW(call eig%setup_matrix(a))

    call a%release(); TEST_IERR()
    call eig%release(); TEST_IERR()

    OUT0("Finished TrilinosEigenSolver_setup_matrix!")

  END_FORTRILINOS_UNIT_TEST(TrilinosEigenSolver_setup_matrix)

  ! -----------------------------setup_matrix_rhs----------------------------- !
    FORTRILINOS_UNIT_TEST(TrilinosEigenSolver_setup_matrix_rhs)
    type(TrilinosEigenSolver) :: eig
    type(TpetraCrsMatrix) :: a
    OUT0("Starting TrilinosEigenSolver_setup_matrix_rhs!")

    call Tpetra_CrsMatrix_CreateIdentity(comm, a)
    eig = TrilinosEigenSolver(comm); TEST_IERR()
    TEST_NOTHROW(call eig%setup_matrix_rhs(a))

    call a%release(); TEST_IERR()
    call eig%release(); TEST_IERR()

    OUT0("Finished TrilinosEigenSolver_setup_matrix_rhs!")

  END_FORTRILINOS_UNIT_TEST(TrilinosEigenSolver_setup_matrix_rhs)

#if 0
  ! ------------------------------setup_operator------------------------------ !
  FORTRILINOS_UNIT_TEST(TrilinosEigenSolver_setup_operator)
    type(TrilinosEigenSolver) :: eig
    type(TpetraMap) :: map
    type(TpetraOperator) :: a
    OUT0("Starting TrilinosEigenSolver_setup_operator!")

    map = TpetraMap(TPETRA_GLOBAL_INVALID, 5, comm)
    a = TpetraOperator(); TEST_IERR()
    !call a%fillComplete()
    eig = TrilinosEigenSolver(comm); TEST_IERR()
    TEST_NOTHROW(call eig%setup_operator(a))

    call a%release(); TEST_IERR()
    call map%release(); TEST_IERR()
    call eig%release(); TEST_IERR()

    OUT0("Finished TrilinosEigenSolver_setup_operator!")

  END_FORTRILINOS_UNIT_TEST(TrilinosEigenSolver_setup_operator)

  ! ----------------------------setup_operator_rhs---------------------------- !
  FORTRILINOS_UNIT_TEST(TrilinosEigenSolver_setup_operator_rhs)
    type(TrilinosEigenSolver) :: eig
    type(TpetraMap) :: map
    type(TpetraOperator) :: a
    OUT0("Starting TrilinosEigenSolver_setup_operator_rhs!")

    map = TpetraMap(TPETRA_GLOBAL_INVALID, 5, comm)
    a = TpetraOperator(); TEST_IERR()
    !call a%fillComplete()
    eig = TrilinosEigenSolver(comm); TEST_IERR()
    TEST_NOTHROW(call eig%setup_operator_rhs(a))

    call a%release(); TEST_IERR()
    call map%release(); TEST_IERR()
    call eig%release(); TEST_IERR()

    OUT0("Finished TrilinosEigenSolver_setup_operator_rhs!")

  END_FORTRILINOS_UNIT_TEST(TrilinosEigenSolver_setup_operator_rhs)
#endif

  ! -------------------------------setup_solver------------------------------- !
  FORTRILINOS_UNIT_TEST(TrilinosEigenSolver_setup_solver)
    type(TrilinosEigenSolver) :: eig
    type(ParameterList) :: plist
    type(TpetraCrsMatrix) :: A
    integer :: num_eigen
    OUT0("Starting TrilinosEigenSolver_setup_solver!")

    plist = ParameterList('Anasazi'); FORTRILINOS_CHECK_IERR()
    TEST_NOTHROW(call load_from_xml(plist, 'davidson.xml'))

    call plist%set('NumEV', num_eigen); FORTRILINOS_CHECK_IERR()
    eig = TrilinosEigenSolver(comm); TEST_IERR()
    call  Tpetra_CrsMatrix_CreateTestMatrix_A(comm,A)
    call A%fillComplete(); FORTRILINOS_CHECK_IERR()
    call eig%setup_matrix(A); TEST_IERR()

    TEST_THROW(call eig%setup_solver(plist))

    call plist%release(); TEST_IERR()
    call A%release(); TEST_IERR()
    call eig%release(); TEST_IERR()

    OUT0("Finished TrilinosEigenSolver_setup_solver!")

  END_FORTRILINOS_UNIT_TEST(TrilinosEigenSolver_setup_solver)
  ! ----------------------------------solve----------------------------------- !
  FORTRILINOS_UNIT_TEST(TrilinosEigenSolver_solve)
    integer :: n, comm_size
    integer :: num_eigen = 1, num_found_eigen, sub_dim, restart_dim
    type(ParameterList) :: plist, dlist
    type(TrilinosEigenSolver) :: eigensolver
    type(TpetraCrsMatrix) :: A
    type(TpetraMultiVector) :: X
    real(scalar_type), dimension(:), allocatable :: evalues
    integer(int_type), dimension(:), allocatable :: eindex

    comm_size = comm%getSize()

    ! Literally copy-paste from fat test, all previous functions required
    ! for this to make any sense
    OUT0("Starting TrilinosEigenSolver_solve!")

    call  Tpetra_CrsMatrix_CreateTestMatrix_A(comm,A)
    call A%fillComplete(); FORTRILINOS_CHECK_IERR()
    n = A%getLocalNumRows()
    num_eigen = 1

    ! This xml file was for another matrix, and we are repurposing but then
    ! rewriting the parameters for our test
    plist = ParameterList('Anasazi'); TEST_IERR()
    call load_from_xml(plist, 'davidson.xml'); TEST_IERR()
    call plist%set("NumEV", num_eigen); FORTRILINOS_CHECK_IERR()
    dlist = plist%sublist('Generalized Davidson')
    sub_dim = (num_eigen+2)*comm%getSize()
    call dlist%set("Maximum Subspace Dimension", sub_dim); FORTRILINOS_CHECK_IERR()
    restart_dim = (num_eigen+1)*comm%getSize()
    call dlist%set("Restart Dimension", restart_dim); FORTRILINOS_CHECK_IERR()

    allocate(evalues(2*num_eigen))
    allocate(eindex(2*num_eigen))

    eigensolver = TrilinosEigenSolver(comm); TEST_IERR()
    call eigensolver%setup_matrix(A); TEST_IERR()
    call eigensolver%setup_solver(plist); TEST_IERR()
    TEST_EQUALITY(num_eigen, eigensolver%max_eigenvalues())
    num_found_eigen = eigensolver%solve(evalues, X, eindex)
    TEST_ASSERT(eigensolver%converged())
    ! This converges in 5 iterations in serial on the docker ci image, 6 on
    ! apple clang with essentially the same installation and options, and 23
    ! iterations when running on 4 processors.
    TEST_ASSERT(eigensolver%num_iters() .gt. 1)
    write(error_unit,*) "converged in iters:", eigensolver%num_iters()

    TEST_EQUALITY(num_found_eigen, num_eigen)
    ! This matrix has lambda=3:
    ! See src/tpetra/test/test_tpetra_crsmatrix_helper.F90
    ! (solver tolerance is 1e-7)
    write(*,*) "evalues:", evalues
    if (abs(3.0_mag_type - evalues(1)) > 1.d-7) then
       write(error_unit,*) "Eigenvalue solver failed"
       stop 1
    endif
    ! Check to make sure the second eigenvalue is 0
    if (abs(0.0_mag_type - evalues(2)) > 1.d-7) then
       write(error_unit,*) "Eigenvalue solver failed"
       stop 1
    endif

    call eigensolver%release(); TEST_IERR()
    call plist%release(); TEST_IERR()
    call dlist%release(); TEST_IERR()
    call X%release(); TEST_IERR()
    call A%release(); TEST_IERR()
    deallocate(eindex, evalues)

    OUT0("Finished TrilinosEigenSolver_solve!")

  END_FORTRILINOS_UNIT_TEST(TrilinosEigenSolver_solve)
end program test_TrilinosEigenSolver
