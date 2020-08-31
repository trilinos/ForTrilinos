! Copyright 2017-2018, UT-Battelle, LLC
!
! SPDX-License-Identifier: BSD-3-Clause
! License-Filename: LICENSE
module myoperators
  use forteuchos
  use fortpetra
  implicit none

  type, extends(ForTpetraOperator) :: TriDiagOperator
    type(TpetraMap) :: row_map, col_map, domain_map, range_map
  contains
    procedure :: apply => my_apply
    procedure :: getDomainMap => my_getDomainMap
    procedure :: getRangeMap => my_getRangeMap
    procedure :: release => delete_TriDiagOperator
  end type
  interface TriDiagOperator
    procedure new_TriDiagOperator
  end interface

contains
  function new_TriDiagOperator(row_map, col_map) &
      result(self)
    use, intrinsic :: ISO_C_BINDING
    type(TpetraMap), intent(in) :: row_map, col_map
    type(TriDiagOperator) :: self
    self = ForTpetraOperator()
    self%row_map = row_map
    self%col_map = col_map
    self%domain_map = row_map
    self%range_map = row_map
  end function

  subroutine my_apply(self, x, y, mode, alpha, beta)
    use, intrinsic :: ISO_C_BINDING
    implicit none
    class(TriDiagOperator), intent(in) :: self
    class(TpetraMultiVector), intent(in) :: x
    class(TpetraMultiVector), intent(in) :: y
    integer(kind(TeuchosETransp)), intent(in) :: mode
    real(scalar_type), intent(in) :: alpha
    real(scalar_type), intent(in) :: beta
    integer :: lid
    integer(global_ordinal_type) :: gid

    type(TeuchosComm) :: comm
    type(TpetraImport) :: import
    type(TpetraMultiVector) :: x_ghosted
    integer :: my_rank, num_procs
    integer :: i, n
    real(scalar_type), dimension(:), pointer :: xdata
    real(scalar_type), dimension(:), pointer :: ydata

    integer, save :: counter = 0

    counter = counter + 1

    comm = self%row_map%getComm()
    my_rank = comm%getRank()
    num_procs = comm%getSize()

    import = TpetraImport(self%domain_map, self%col_map)
    x_ghosted = TpetraMultiVector(self%col_map, 1_size_type)
    call x_ghosted%doImport(x, import, TpetraINSERT)
    call import%release()

    xdata => x_ghosted%getData        (1_size_type)
    ydata => y        %getDataNonConst(1_size_type)
    n = y%getLocalLength()

    ! Sometimes, ydata may be unitialized (when beta is 0), potentially containing
    ! signaling NaNs. Therefore, for beta = 0, we explicitly zero it out.
    if (beta .eq. 0) then
      do i = 1, n
        ydata(i) = 0
      end do
    else
      do i = 1, n
        ydata(i) = beta * ydata(i)
      end do
    end if

    ! y = alpha * A*x + beta * y
    do i = 1, n
      gid = self%range_map%getGlobalElement(i)

      ! A has [-1 2 -1] stencil
      if (i > 1 .or. my_rank > 0) then
        lid = self%col_map%getLocalElement(gid-1)
        ydata(i) = ydata(i) - alpha*xdata(lid)
      end if

      lid = self%col_map%getLocalElement(gid)
      ydata(i) = ydata(i) + 2*alpha*xdata(lid)

      if (i < n .or. my_rank .ne. num_procs-1) then
        lid = self%col_map%getLocalElement(gid+1)
        ydata(i) = ydata(i) - alpha*xdata(lid)
      end if
    end do

    nullify(xdata)
    nullify(ydata)

    call x_ghosted%release()
    call comm%release()

  end subroutine

  function my_getDomainMap(self) &
      result(domain_map)
    use, intrinsic :: ISO_C_BINDING
    implicit none
    class(TriDiagOperator), intent(in) :: self
    type(TpetraMap) :: domain_map

    domain_map = self%domain_map
  end function

  function my_getRangeMap(self) &
      result(range_map)
    use, intrinsic :: ISO_C_BINDING
    implicit none
    class(TriDiagOperator), intent(in) :: self
    type(TpetraMap) :: range_map

    range_map = self%range_map
  end function

  subroutine delete_TriDiagOperator(self)
    class(TriDiagOperator), intent(inout) :: self

    call self%row_map%release()
    call self%col_map%release()
    call self%domain_map%release()
    call self%range_map%release()

#ifdef __GNUC__
    ! FIXME This segfaults with Flang
    ! Call base class release()
    call self%ForTpetraOperator%release()
#endif
  end subroutine


end module


program main

#include "ForTrilinos_config.h"
#include "ForTrilinos.h"

  use ISO_FORTRAN_ENV
  use, intrinsic :: ISO_C_BINDING

#if FORTRILINOS_USE_MPI
  use mpi
#endif

#include "ForTrilinos.h"
  use fortrilinos_hl
  use forteuchos
  use fortpetra
  use fortest
  use myoperators
  implicit none

  integer :: my_rank, num_procs

  integer(global_size_type) :: n_global
  integer(size_type) :: max_entries_per_row, num_vecs = 1, lda
  integer(int_type) :: row_nnz

  integer :: ierr
  integer :: i, n
  integer(global_ordinal_type) :: offset

  type(TeuchosComm) :: comm
  type(ParameterList) :: plist, linear_solver_list, belos_list, solver_list, krylov_list
  type(TrilinosSolver) :: solver_handle
  type(TpetraMap) :: map
  type(TpetraCrsMatrix) :: A
  type(TpetraMultiVector) :: B, X, residual
  class(ForTpetraOperator), allocatable :: op

  real(scalar_type), dimension(:), allocatable :: lhs, rhs
  real(norm_type), dimension(:), allocatable :: norms
  integer(global_ordinal_type), dimension(:), allocatable :: cols
  real(scalar_type), dimension(:), allocatable :: vals
  real(scalar_type) :: r0, sone = 1., szero = 0., tol

  n = 10

#if FORTRILINOS_USE_MPI
  ! Initialize MPI subsystem
  call MPI_INIT(ierr)
  if (ierr /= 0) then
    write(*,*) "MPI failed to init"
    stop 1
  endif

  comm = TeuchosComm(MPI_COMM_WORLD)
#else
  comm = TeuchosComm()
#endif

  my_rank = comm%getRank()
  num_procs = comm%getSize()

  ! Read in the parameterList
  plist = ParameterList("Stratimikos"); FORTRILINOS_CHECK_IERR()
  call load_from_xml(plist, "stratimikos.xml"); FORTRILINOS_CHECK_IERR()

  ! Get tolerance from the parameter list
  linear_solver_list = plist%sublist('Linear Solver Types')
  belos_list = linear_solver_list%sublist(plist%get_string('Linear Solver Type'))
  solver_list = belos_list%sublist('Solver Types')
  krylov_list = solver_list%sublist(belos_list%get_string('Solver Type'))

  tol = krylov_list%get_real("Convergence Tolerance")

  ! ------------------------------------------------------------------
  ! Step 0: Construct tri-diagonal matrix
  n_global = -1
  map = TpetraMap(n_global, n, comm); FORTRILINOS_CHECK_IERR()

  max_entries_per_row = 3
  A = TpetraCrsMatrix(map, max_entries_per_row, TpetraStaticProfile)

  allocate(cols(max_entries_per_row))
  allocate(vals(max_entries_per_row))
  offset = n * my_rank
  do i = 1, n
    row_nnz = 1
    if (i .ne. 1 .or. my_rank > 0) then
      cols(row_nnz) = offset + i-1
      vals(row_nnz) = -1.0
      row_nnz = row_nnz + 1
    end if
    cols(row_nnz) = offset + i
    vals(row_nnz) = 2.0
    row_nnz = row_nnz + 1
    if (i .ne. n .or. my_rank .ne. num_procs-1) then
      cols(row_nnz) = offset + i+1
      vals(row_nnz) = -1.0
      row_nnz = row_nnz + 1
    end if

    call A%insertGlobalValues(offset + i, cols(1:row_nnz-1), vals(1:row_nnz-1)); FORTRILINOS_CHECK_IERR()
  end do
  call A%fillComplete(); FORTRILINOS_CHECK_IERR()

  ! The solution X(i) = i-1
  allocate(lhs(n))
  allocate(rhs(n))
  if (my_rank > 0) then
    rhs(1) = 0.0
  else
    rhs(1) = -1.0
  end if
  if (my_rank .ne. num_procs-1) then
    rhs(n) = 0.0
  else
    rhs(n) = offset+n
  end if
  do i = 2, n-1
    rhs(i) = 0.0
  end do
  lda = n

  B = TpetraMultiVector(map, rhs, lda, num_vecs); FORTRILINOS_CHECK_IERR()
  X = TpetraMultiVector(map, num_vecs); FORTRILINOS_CHECK_IERR()
  residual = TpetraMultiVector(map, num_vecs, .false.); FORTRILINOS_CHECK_IERR()

  allocate(norms(1))

  ! Step 0: create a handle
  solver_handle = TrilinosSolver(); FORTRILINOS_CHECK_IERR()

  ! ------------------------------------------------------------------
  ! Explicit setup and solve
  ! ------------------------------------------------------------------

  ! Use case 1: setup and solve the problem

  ! Step 1: initialize a handle
  call solver_handle%init(comm); FORTRILINOS_CHECK_IERR()

  ! Step 2: setup the problem
  call solver_handle%setup_matrix(A); FORTRILINOS_CHECK_IERR()

  ! Step 3: setup the solver
  call solver_handle%setup_solver(plist); FORTRILINOS_CHECK_IERR()

  ! Step 4: solve the system
  call X%randomize()
  ! Calculate initial residual
  call A%apply(X, residual, TeuchosNO_TRANS, sone, szero); FORTRILINOS_CHECK_IERR()
  call residual%update(sone, B, -sone); FORTRILINOS_CHECK_IERR()
  call residual%norm2(norms); FORTRILINOS_CHECK_IERR()
  r0 = norms(1)
  call solver_handle%solve(B, X); FORTRILINOS_CHECK_IERR()

  ! Check the solution
  call A%apply(X, residual, TeuchosNO_TRANS, sone, szero); FORTRILINOS_CHECK_IERR()
  call residual%update(sone, B, -sone); FORTRILINOS_CHECK_IERR()
  call residual%norm2(norms); FORTRILINOS_CHECK_IERR()

  if (norms(1)/r0 > tol) then
    write(error_unit, '(A)') 'The solver did not converge to the specified residual!'
    stop 1
  end if

  ! Use case 2: Reuse preconditioner with a different matrix
  call A%resumeFill()
  do i = 1, n
    cols(1) = offset + i
    vals(1) = 3.0
    row_nnz = A%replaceGlobalValues(offset + i, cols(1:1), vals(1:1)); FORTRILINOS_CHECK_IERR()
  end do
  call A%fillComplete(); FORTRILINOS_CHECK_IERR()

  ! Step 2a: setup the problem that has already been setup
  call solver_handle%setup_matrix(A); FORTRILINOS_CHECK_IERR()

  ! Step 4a: solve the system with the original preconditioner
  call X%randomize()
  ! Calculate initial residual
  call A%apply(X, residual, TeuchosNO_TRANS, sone, szero); FORTRILINOS_CHECK_IERR()
  call residual%update(sone, B, -sone); FORTRILINOS_CHECK_IERR()
  call residual%norm2(norms); FORTRILINOS_CHECK_IERR()
  r0 = norms(1)
  call solver_handle%solve(B, X); FORTRILINOS_CHECK_IERR()

  ! Check the solution
  call A%apply(X, residual, TeuchosNO_TRANS, sone, szero); FORTRILINOS_CHECK_IERR()
  call residual%update(sone, B, -sone); FORTRILINOS_CHECK_IERR()
  call residual%norm2(norms); FORTRILINOS_CHECK_IERR()

  if (norms(1)/r0 > tol) then
    write(error_unit, '(A)') 'The solver did not converge to the specified residual!'
    stop 1
  end if

  ! Step 5: clean up
  call solver_handle%finalize(); FORTRILINOS_CHECK_IERR()


  ! Change values back for implicit comparison
  call A%resumeFill()
  do i = 1, n
    cols(1) = offset + i
    vals(1) = 2.0
    row_nnz = A%replaceGlobalValues(offset + i, cols(1:1), vals(1:1)); FORTRILINOS_CHECK_IERR()
  end do
  call A%fillComplete(); FORTRILINOS_CHECK_IERR()

  ! ------------------------------------------------------------------
  ! Implicit (inversion-of-control) setup [ no solve ]
  ! ------------------------------------------------------------------
  ! We cannot use most preconditioners without a matrix, so we remove any from
  ! the parameter list. We also adjust the number of iterations so that it is
  ! sufficient for convergence
  call plist%set('Preconditioner Type', 'None')
  call krylov_list%set('Maximum Iterations', 333)

  allocate(op, source=TriDiagOperator(map, A%getColMap()))
  call init_ForTpetraOperator(op); FORTRILINOS_CHECK_IERR()

  ! Step 1: initialize a handle
  call solver_handle%init(comm); FORTRILINOS_CHECK_IERR()

  ! Step 2: setup the problem
  ! Implicit (inversion-of-control) setup
  call solver_handle%setup_operator(op); FORTRILINOS_CHECK_IERR()

  ! Step 3: setup the solver

  call solver_handle%setup_solver(plist); FORTRILINOS_CHECK_IERR()

  call krylov_list%release()

  ! Step 4: solve the system
  call X%randomize()
  ! Calculate initial residual
  call A%apply(X, residual, TeuchosNO_TRANS, sone, szero); FORTRILINOS_CHECK_IERR()
  call residual%update(sone, B, -sone); FORTRILINOS_CHECK_IERR()
  call residual%norm2(norms); FORTRILINOS_CHECK_IERR()
  r0 = norms(1)
  call solver_handle%solve(B, X); FORTRILINOS_CHECK_IERR()

  ! Check the solution
  call A%apply(X, residual, TeuchosNO_TRANS, sone, szero); FORTRILINOS_CHECK_IERR()
  call residual%update(sone, B, -sone); FORTRILINOS_CHECK_IERR()
  call residual%norm2(norms); FORTRILINOS_CHECK_IERR()
  if (norms(1)/r0 > tol) then
    write(error_unit, '(A)') 'The solver did not converge to the specified residual!'
    stop 666
  end if

  ! Step 5: clean up
  call solver_handle%finalize(); FORTRILINOS_CHECK_IERR()

  call krylov_list%release; FORTRILINOS_CHECK_IERR()
  call solver_list%release; FORTRILINOS_CHECK_IERR()
  call belos_list%release; FORTRILINOS_CHECK_IERR()
  call linear_solver_list%release; FORTRILINOS_CHECK_IERR()
  call op%release(); FORTRILINOS_CHECK_IERR()
  deallocate(op)
  ! ------------------------------------------------------------------

  call solver_handle%release(); FORTRILINOS_CHECK_IERR()
  call plist%release(); FORTRILINOS_CHECK_IERR()
  call X%release(); FORTRILINOS_CHECK_IERR()
  call B%release(); FORTRILINOS_CHECK_IERR()
  call A%release(); FORTRILINOS_CHECK_IERR()
  call map%release(); FORTRILINOS_CHECK_IERR()
  call comm%release(); FORTRILINOS_CHECK_IERR()
  deallocate(norms)
  deallocate(cols)
  deallocate(vals)
  deallocate(lhs)
  deallocate(rhs)

#if FORTRILINOS_USE_MPI
  ! Finalize MPI must be called after releasing all handles
  call MPI_FINALIZE(ierr)
#endif


end program
