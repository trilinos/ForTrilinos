! Copyright 2017, UT-Battelle, LLC
!
! SPDX-License-Identifier: BSD-3-Clause
! License-Filename: LICENSE
module x_eigen_handle
  use iso_c_binding
  implicit none
contains
  subroutine matvec(n, x, y) BIND(C)
    use, intrinsic :: ISO_C_BINDING
    integer(c_int), intent(in), value :: n
    real(c_double), dimension(:), intent(in) :: x(*)
    real(c_double), dimension(:), intent(out) :: y(*)

    integer(c_int) :: i

    ! dummy operator
    do i = 1, n
      y(i) = x(i)
    end do

  end subroutine
end module x_eigen_handle

program main

#include "FortranTestMacros.h"
#include "ForTrilinosSimpleInterface_config.hpp"


  use ISO_FORTRAN_ENV
  use, intrinsic :: ISO_C_BINDING
  use fortrilinos
  use x_eigen_handle
  use forteuchos
#ifdef HAVE_MPI
  use mpi
#endif
  implicit none

  integer(int_type) :: i
  integer(int_type) :: n, nnz;
  integer(int_type) :: my_rank, num_procs

  integer(local_ordinal_type), dimension(:), allocatable :: row_ptrs
  integer(global_ordinal_type), dimension(:), allocatable :: row_inds, col_inds
  real(scalar_type), dimension(:), allocatable :: values, evalues, evectors

  integer(global_ordinal_type) :: cur_pos, offset
  real(norm_type) :: norm

  type(ParameterList) :: plist
  type(EigenHandle) :: tri_handle
  type(TeuchosComm) :: comm

  n = 50
  nnz = 3*n

  my_rank = 0
  num_procs = 1

#ifdef HAVE_MPI
  ! Initialize MPI subsystem
  call MPI_INIT(ierr)
  if (ierr /= 0) then
    write(*,*) "MPI failed to init"
    stop 1
  endif

  call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
  EXPECT_EQ(0, ierr)

  call comm%create(MPI_COMM_WORLD)
#else
  call comm%create()
#endif

  ! Read in the parameterList
  call plist%create("Anasazi")
  call load_from_xml(plist, "davidson.xml")

  if (ierr /= 0) then
    write(*,*) "Got error ", ierr, ":", trim(serr)
    stop 1
  endif

  ! ------------------------------------------------------------------
  ! Step 0: Construct tri-diagonal matrix
  allocate(row_inds(n))
  allocate(row_ptrs(n+1))
  allocate(col_inds(nnz))
  allocate(values(nnz))
  row_ptrs(1) = 0
  cur_pos = 1
  offset  = n * my_rank
  do i = 1, n
    ! Is logic evaluated the same way in fortran?
    if (i .ne. 1 .or. my_rank > 0) then
      col_inds(cur_pos) = offset + i-1
      values  (cur_pos) = -1.0
      cur_pos = cur_pos + 1
    end if
    col_inds(cur_pos) = offset + i
    values  (cur_pos) = 2.0
    cur_pos = cur_pos + 1
    if (i .ne. n .or. my_rank .ne. num_procs-1) then
      col_inds(cur_pos) = offset + i+1
      values  (cur_pos) = -1.0
      cur_pos = cur_pos + 1
    end if
    row_ptrs(i+1) = cur_pos-1;

    row_inds(i) = offset + i
  end do
  nnz = cur_pos-1

  ! The solution
  allocate(evalues(1))
  allocate(evectors(n))

  ! Step 0.5: crate a handle
  call tri_handle%create()
  if (ierr /= 0) then
    write(*,*) "Got error ", ierr, ": ", trim(serr)
    stop 1
  endif

  ! ------------------------------------------------------------------
  ! Explicit setup and solve
  ! ------------------------------------------------------------------
  ! Step 1: initialize a handle
  call tri_handle%init(comm)
  if (ierr /= 0) then
    write(*,*) "Got error ", ierr, ": ", trim(serr)
    stop 1
  endif

  ! Step 2: setup the problem
  call tri_handle%setup_matrix(row_inds, row_ptrs, col_inds, values)
  if (ierr /= 0) then
    write(*,*) "Got error ", ierr, ": ", trim(serr)
    stop 1
  endif

  ! Step 3: setup the solver
  call tri_handle%setup_solver(plist)
  if (ierr /= 0) then
    write(*,*) "Got error ", ierr, ": ", trim(serr)
    stop 1
  endif

  ! Step 4: solve the system
  call tri_handle%solve(evalues, evectors)
  if (ierr /= 0) then
    write(*,*) "Got error ", ierr, ": ", trim(serr)
    stop 1
  endif

  ! Check the solution


  ! Step 5: clean up
  call tri_handle%finalize()
  if (ierr /= 0) then
    write(*,*) "Got error ", ierr, ": ", trim(serr)
    stop 1
  endif

#if 0
  ! ------------------------------------------------------------------
  ! Implicit (inversion-of-control) setup [ no solve ]
  ! ------------------------------------------------------------------
  ! Step 1: initialize a handle
  call tri_handle%init(comm)
  if (ierr /= 0) then
    write(*,*) "Got error ", ierr, ": ", trim(serr)
    stop 1
  endif

  ! Step 2: setup the problem
  ! Implicit (inversion-of-control) setup
  call tri_handle%setup_operator(row_inds, C_FUNLOC(matvec))
  if (ierr /= 0) then
    write(*,*) "Got error ", ierr, ": ", trim(serr)
    stop 1
  endif

  ! Step 3: setup the solver
  ! We cannot use most preconditioners without a matrix, so
  ! we remove any from the parameter list
  call plist%set("Preconditioner Type", "None")
  call tri_handle%setup_solver(plist)
  if (ierr /= 0) then
    write(*,*) "Got error ", ierr, ": ", trim(serr)
    stop 1
  endif

  ! Step 4: solve the system
  ! We only check that it runs, but do not check the result as
  ! we are using a dummy operator
  call tri_handle%solve(rhs, lhs)
  if (ierr /= 0) then
    write(*,*) "Got error ", ierr, ":", trim(serr)
    stop 1
  endif

  ! Step 5: clean up
  call tri_handle%finalize()
  if (ierr /= 0) then
    write(*,*) "Got error ", ierr, ":", trim(serr)
    stop 1
  endif

  ! ------------------------------------------------------------------

  call plist%release()
  call tri_handle%release()
  if (ierr /= 0) then
    write(*,*) "Got error ", ierr, ":", trim(serr)
    stop 1
  endif
#endif

  call comm%release()

#ifdef HAVE_MPI
  ! Finalize MPI must be called after releasing all handles
  call MPI_FINALIZE(ierr)
  EXPECT_EQ(0, ierr)
#endif


end program
