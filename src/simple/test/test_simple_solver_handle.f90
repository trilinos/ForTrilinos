! Copyright 2017-2018, UT-Battelle, LLC
!
! SPDX-License-Identifier: BSD-3-Clause
! License-Filename: LICENSE
module myoperators
  use forteuchos
  use fortpetra
  implicit none

  type, extends(ForTpetraOperator) :: TriDiagOperator
    type(TpetraMap) :: domain_map
    type(TpetraMap) :: range_map
  contains
    procedure :: apply => my_apply
    procedure :: getDomainMap => my_getDomainMap
    procedure :: getRangeMap => my_getRangeMap
  end type
  interface TriDiagOperator
    procedure new_TriDiagOperator
  end interface

contains
  function new_TriDiagOperator(domain_map, range_map) &
      result(self)
    use, intrinsic :: ISO_C_BINDING
    type(TpetraMap), intent(in) :: domain_map
    type(TpetraMap), intent(in) :: range_map
    type(TriDiagOperator) :: self
    write(*,*) "Called constructor"
    self%domain_map = domain_map
    self%range_map = range_map
  end function

  subroutine my_apply(self, x, y, mode, alpha, beta)
    use, intrinsic :: ISO_C_BINDING
    implicit none
    class(TriDiagOperator), intent(in) :: self
    class(TpetraMultiVector), intent(in) :: x
    class(TpetraMultiVector), intent(inout) :: y
    integer(kind(TeuchosETransp)), intent(in) :: mode
    real(scalar_type), intent(in) :: alpha
    real(scalar_type), intent(in) :: beta
    integer(size_type) :: n, i
    real(scalar_type), dimension(:), pointer :: xdata
    real(scalar_type), dimension(:), pointer :: ydata

    xdata => x%getData(INT(1,size_type))
    ydata => y%getData(INT(1,size_type))
    n = x%getLocalLength()

    ! y = alpha * A*x + beta * y
    do i = 1, n
      ydata(i) = beta * ydata(i)

      ! A has [-1 2 -1] stencil
      ydata(i) = ydata(i) + 2*alpha*xdata(i)
      if (i > 1) then
        ydata(i) = ydata(i) - alpha*xdata(i-1)
      end if
      if (i < n) then
        ydata(i) = ydata(i) - alpha*xdata(i+1)
      end if
    end do

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

end module


program main

#include "ForTrilinosSimpleInterface_config.hpp"
#include "ForTrilinos.h"

  use ISO_FORTRAN_ENV
  use, intrinsic :: ISO_C_BINDING

#ifdef HAVE_MPI
  use mpi
#endif

#include "ForTrilinos.h"
  use fortrilinos
  use forteuchos
  use fortpetra
  use fortest
  use myoperators
  implicit none

  integer(int_type) :: my_rank, num_procs

  integer(global_size_type) :: n_global
  integer(size_type) :: n, max_entries_per_row, num_vecs = 1, lda
  integer(int_type) :: row_nnz

  integer :: ierr
  integer(local_ordinal_type) :: i
  integer(global_ordinal_type) :: offset
  real(scalar_type) :: one = 1.0
  real(norm_type) :: norm

  type(TeuchosComm) :: comm
  type(ParameterList) :: plist
  type(TrilinosSolver) :: solver_handle
  type(TpetraMap) :: map
  type(TpetraCrsMatrix) :: A
  type(TpetraMultiVector) :: B, X, Xtrue
  class(ForTpetraOperator), allocatable :: op

  real(scalar_type), dimension(:), allocatable :: lhs, rhs
  real(norm_type), dimension(:), allocatable :: norms
  integer(global_ordinal_type), dimension(:), allocatable :: cols
  real(scalar_type), dimension(:), allocatable :: vals

  n = 50

#ifdef HAVE_MPI
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

  write(*,*) "Processor ", my_rank, " of ", num_procs

  ! Read in the parameterList
  plist = ParameterList("Stratimikos"); FORTRILINOS_CHECK_IERR()
  call load_from_xml(plist, "stratimikos.xml"); FORTRILINOS_CHECK_IERR()

  ! ------------------------------------------------------------------
  ! Step 0: Construct tri-diagonal matrix
  n_global = -1
  map = TpetraMap(n_global, n, comm); FORTRILINOS_CHECK_IERR()

  max_entries_per_row = 3
  A = TpetraCrsMatrix(map, max_entries_per_row, TpetraDynamicProfile)

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

  ! This automatically zeroes out X
  X = TpetraMultiVector(map, num_vecs); FORTRILINOS_CHECK_IERR()

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
  do i = 1, n
    lhs(i) = offset + i-1
  end do
  lda = n

  Xtrue = TpetraMultiVector(map, lhs, lda, num_vecs); FORTRILINOS_CHECK_IERR()
  B = TpetraMultiVector(map, rhs, lda, num_vecs); FORTRILINOS_CHECK_IERR()

  ! Step 0: create a handle
  solver_handle = TrilinosSolver(); FORTRILINOS_CHECK_IERR()

  ! ------------------------------------------------------------------
  ! Explicit setup and solve
  ! ------------------------------------------------------------------

  ! Step 1: initialize a handle
  call solver_handle%init(comm); FORTRILINOS_CHECK_IERR()

  ! Step 2: setup the problem
  call solver_handle%setup_matrix(A); FORTRILINOS_CHECK_IERR()

  ! Step 3: setup the solver
  call plist%set("Preconditioner Type", "None")
  call solver_handle%setup_solver(plist); FORTRILINOS_CHECK_IERR()

  ! Step 4: solve the system
  call X%randomize()
  call solver_handle%solve(B, X); FORTRILINOS_CHECK_IERR()

  ! Check the solution
  allocate(norms(1))
  call X%update(-one, Xtrue, one); FORTRILINOS_CHECK_IERR()
  call X%norm2(norms); FORTRILINOS_CHECK_IERR()

  ! TODO: Get the tolerance out of the parameter list
  if (norms(1) > 1e-6) then
    write(error_unit, '(A)') 'The solver did not converge to the specified residual!'
    stop 1
  end if

  ! Step 5: clean up
  call solver_handle%finalize(); FORTRILINOS_CHECK_IERR()

  ! ------------------------------------------------------------------
  ! Implicit (inversion-of-control) setup [ no solve ]
  ! ------------------------------------------------------------------
  allocate(op, source=TriDiagOperator(map, map))
  ! Step 1: initialize a handle
  call solver_handle%init(comm); FORTRILINOS_CHECK_IERR()

  ! Step 2: setup the problem
  ! Implicit (inversion-of-control) setup
  call init_ForTpetraOperator(op); FORTRILINOS_CHECK_IERR()
  call solver_handle%setup_operator(op); FORTRILINOS_CHECK_IERR()

  ! Step 3: setup the solver
  ! We cannot use most preconditioners without a matrix, so
  ! we remove any from the parameter list
  call plist%set("Preconditioner Type", "None")
  call solver_handle%setup_solver(plist); FORTRILINOS_CHECK_IERR()

  ! Step 4: solve the system
  call X%randomize()
  call solver_handle%solve(B, X); FORTRILINOS_CHECK_IERR()

  call X%update(-one, Xtrue, one); FORTRILINOS_CHECK_IERR()
  call X%norm2(norms); FORTRILINOS_CHECK_IERR()
  if (norms(1) > 1e-10) then
    write(error_unit, '(A)') 'The implicit result differs from explicit!'
    write(error_unit, '(A, ES14.7)') '   diff_norm =', norms(1)
    stop 666
  end if

  ! Step 5: clean up
  call solver_handle%finalize(); FORTRILINOS_CHECK_IERR()

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

#ifdef HAVE_MPI
  ! Finalize MPI must be called after releasing all handles
  call MPI_FINALIZE(ierr)
#endif


end program
