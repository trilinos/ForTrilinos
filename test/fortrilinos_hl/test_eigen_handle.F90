! Copyright 2017-2018, UT-Battelle, LLC
!
! SPDX-License-Identifier: BSD-3-Clause
! License-Filename: LICENSE
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

  implicit none

  integer :: my_rank, num_procs
  integer(global_size_type) :: n_global
  integer:: i, n
  integer(size_type) :: max_entries_per_row
  integer :: num_eigen = 1, num_found_eigen
  integer(int_type) :: row_nnz

  integer :: ierr
  integer(global_ordinal_type) :: offset

  type(TeuchosComm) :: comm
  type(ParameterList) :: plist
  type(TrilinosEigenSolver) :: eigen_handle
  type(TpetraMap) :: map
  type(TpetraCrsMatrix) :: A
  type(TpetraMultiVector) :: X

  real(scalar_type), dimension(:), allocatable :: evalues
  integer(int_type), dimension(:), allocatable :: eindex
  integer(global_ordinal_type), dimension(:), allocatable :: cols
  real(scalar_type), dimension(:), allocatable :: vals

  n = 50

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

  write(*,*) "Processor ", my_rank, " of ", num_procs

  ! Read in the parameterList
  plist = ParameterList("Anasazi"); FORTRILINOS_CHECK_IERR()
  call load_from_xml(plist, "davidson.xml"); FORTRILINOS_CHECK_IERR()

  call plist%set("NumEV", num_eigen); FORTRILINOS_CHECK_IERR()

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

  ! The solution
  allocate(evalues(2*num_eigen))
  allocate(eindex(2*num_eigen))

  ! Step 0: crate a handle
  eigen_handle = TrilinosEigenSolver(); FORTRILINOS_CHECK_IERR()

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
  num_found_eigen = eigen_handle%solve(evalues, X, eindex); FORTRILINOS_CHECK_IERR()

  ! FIXME: Check the solution
  if (num_found_eigen < num_eigen) then
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
  deallocate(eindex)
  deallocate(evalues)
  deallocate(cols)
  deallocate(vals)

#if FORTRILINOS_USE_MPI
  ! Finalize MPI must be called after releasing all handles
  call MPI_FINALIZE(ierr)
#endif


end program
