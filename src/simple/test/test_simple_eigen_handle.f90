! Copyright 2017, UT-Battelle, LLC
!
! SPDX-License-Identifier: BSD-3-Clause
! License-Filename: LICENSE
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

  implicit none

  integer(int_type) :: my_rank, num_procs

  integer(global_size_type) :: n_global
  integer(size_type) :: n, max_entries_per_row, lda, num_eigen = 1
  integer(int_type) :: num_eigen_int
  integer(int_type) :: row_nnz

  integer :: ierr
  integer(local_ordinal_type) :: i
  integer(global_ordinal_type) :: offset
  real(scalar_type) :: one = 1.0

  type(TeuchosComm) :: comm
  type(ParameterList) :: plist
  type(TrilinosEigenSolver) :: eigen_handle
  type(TpetraMap) :: map
  type(TpetraCrsMatrix) :: A
  type(TpetraMultiVector) :: X

  real(scalar_type), dimension(:), allocatable :: evalues, evectors
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

  call comm%create(MPI_COMM_WORLD)
#else
  call comm%create()
#endif

  my_rank = comm%getRank()
  num_procs = comm%getSize()

  write(*,*) "Processor ", my_rank, " of ", num_procs

  ! Read in the parameterList
  call plist%create("Anasazi"); FORTRILINOS_CHECK_IERR()
  call load_from_xml(plist, "davidson.xml"); FORTRILINOS_CHECK_IERR()

  num_eigen_int = num_eigen
  call plist%set("NumEV", num_eigen_int); FORTRILINOS_CHECK_IERR()

  ! ------------------------------------------------------------------
  ! Step 0: Construct tri-diagonal matrix
  n_global = -1
  call map%create(n_global, n, comm); FORTRILINOS_CHECK_IERR()

  max_entries_per_row = 3
  call A%create(map, max_entries_per_row, TpetraDynamicProfile)

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

  ! The solution
  allocate(evalues(num_eigen))
  call X%create(map, num_eigen)

  ! Step 0: crate a handle
  call eigen_handle%create(); FORTRILINOS_CHECK_IERR()

  ! ------------------------------------------------------------------
  ! Explicit setup and solve
  ! ------------------------------------------------------------------
  ! Step 1: initialize a handle
  call eigen_handle%init(comm); FORTRILINOS_CHECK_IERR()

  ! Step 2: setup the problem
  call eigen_handle%setup_matrix(A); FORTRILINOS_CHECK_IERR()

  ! Step 3: setup the solver
  call eigen_handle%setup_solver(plist); FORTRILINOS_CHECK_IERR()

  ! Step 4: solve the system
  call eigen_handle%solve(evalues, X); FORTRILINOS_CHECK_IERR()

  ! FIXME: Check the solution
  if (size(evalues) /= num_eigen) then
    write(error_unit, '(A)') 'The number of returned eigenvalues does not match!'
    stop 1
  end if
  write(*,*) "Computed eigenvalues: ", evalues(1)

  ! Step 5: clean up
  call eigen_handle%finalize(); FORTRILINOS_CHECK_IERR()

  call eigen_handle%release(); FORTRILINOS_CHECK_IERR()
  call plist%release(); FORTRILINOS_CHECK_IERR()
  call X%release(); FORTRILINOS_CHECK_IERR()
  call A%release(); FORTRILINOS_CHECK_IERR()
  call map%release(); FORTRILINOS_CHECK_IERR()
  call comm%release(); FORTRILINOS_CHECK_IERR()
  deallocate(evalues)
  deallocate(cols)
  deallocate(vals)

#ifdef HAVE_MPI
  ! Finalize MPI must be called after releasing all handles
  call MPI_FINALIZE(ierr)
#endif


end program
