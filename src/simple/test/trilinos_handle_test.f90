program main

#include "FortranTestMacros.h"

  use ISO_FORTRAN_ENV
  use, intrinsic :: ISO_C_BINDING
  use simpleinterface
  use forteuchos
  use mpi
  implicit none

  integer(c_int) :: i
  integer(c_int) :: n, nnz;
  integer(c_int) :: my_rank, num_procs

  integer(c_int), dimension(:), allocatable :: row_inds, col_inds, row_ptrs
  real(c_double), dimension(:), allocatable :: lhs, rhs, values

  integer(c_int) :: cur_pos, offset
  real(c_double) :: norm

  type(ParameterList) :: plist
  type(TrilinosHandle) :: tri_handle

  integer :: ierr

  n = 50
  nnz = 3*n

  my_rank = 0
  num_procs = 1

  ! Initialize MPI subsystem
  call MPI_INIT(ierr)
  EXPECT_EQ(ierr, 0)

  call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank,   ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
  EXPECT_EQ(ierr, 0)

  ! Read in the parameterList
  call plist%create("Stratimikos")
  call load_from_xml(plist, "stratimikos.xml")

  ! Step 0: Construct tri-diagonal matrix, and rhs
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

  allocate(lhs(n))
  do i = 1, n
    lhs(i) = 0.0
  end do

  ! The solution lhs(i) = i-1
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

  ! Step 0: crate a handle
  call tri_handle%create()

  ! Step 1: initialize a handle
  ! call init(ierr)
  call tri_handle%init(MPI_COMM_WORLD, ierr)
  EXPECT_EQ(ierr, 0)

  ! Step 2: setup the problem
  call tri_handle%setup_matrix(n, row_inds, row_ptrs, nnz, col_inds, values, ierr)
  EXPECT_EQ(ierr, 0)

  ! // Step 3: setup the solver
  call tri_handle%setup_solver(plist, ierr)
  EXPECT_EQ(ierr, 0)

  ! // Step 4: solve the system
  call tri_handle%solve(n, rhs, lhs, ierr)
  EXPECT_EQ(ierr, 0)

  ! Check the solution
  norm = 0.0
  do i = 1, n
    norm = norm + (lhs(i) - (offset+i-1))*(lhs(i) - (offset+i-1));
  end do
  norm = sqrt(norm);

  write(0,*) "norm = ", norm

  ! TODO: Get the tolerance out of the parameter list
  EXPECT_TRUE(norm < 1e-6)


  ! Step 5: clean up
  call tri_handle%finalize(ierr)
  EXPECT_EQ(ierr, 0)

  call MPI_FINALIZE(ierr)
  EXPECT_EQ(ierr, 0)

  call plist%release()
  call tri_handle%release()

end program
