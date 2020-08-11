
! Copyright 2017, UT-Battelle, LLC
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

  integer(int_type) :: my_rank, num_procs

  integer(global_size_type) :: n_global
  integer(size_type) :: max_entries_per_row
  integer :: n
  integer(size_type) :: num_eigen = 1, num_found_eigen

  integer :: ierr

  type(TeuchosComm) :: comm
  type(ParameterList) :: plist
  type(TrilinosEigenSolver) :: eigen_handle
  type(TpetraMap) :: map
  type(TpetraCrsMatrix) :: A,B
  type(TpetraMultiVector) :: X

  real(scalar_type), dimension(:), allocatable :: evalues
  integer(int_type), dimension(:), allocatable :: eindex
  integer(global_ordinal_type), dimension(:) :: cols(2)
  real(scalar_type), dimension(:) :: vals(2)

  n = 30

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
  call plist%set("Preconditioner Type", "IDENTITY")

  call plist%set("NumEV", 1); FORTRILINOS_CHECK_IERR()

  ! ------------------------------------------------------------------
  ! Step 0: Construct A and B matrices
  n_global = -1
  map = TpetraMap(n_global, n, comm); FORTRILINOS_CHECK_IERR()

  max_entries_per_row = 2
  A = TpetraCrsMatrix(map, map, max_entries_per_row)
  cols(1) = 1
  cols(2) = 2
  vals(1) = 0.0015
  vals(2) = 0.325
  call A%insertGlobalValues(1_global_ordinal_type, cols(1:2), vals(1:2)); FORTRILINOS_CHECK_IERR()

  B = TpetraCrsMatrix(map, map, max_entries_per_row)
  cols(1) = 1
  vals(1) = 0.1208
  call B%insertGlobalValues(1_global_ordinal_type, cols(1:1), vals(1:1)); FORTRILINOS_CHECK_IERR()
  cols(1) = 1
  cols(2) = 2
  vals(1) = -0.117
  vals(2) = 0.184
  call B%insertGlobalValues(2_global_ordinal_type, cols(1:2), vals(1:2)); FORTRILINOS_CHECK_IERR()

  call A%fillComplete(); FORTRILINOS_CHECK_IERR()
  call B%fillComplete(); FORTRILINOS_CHECK_IERR()

  ! The solution
  allocate(evalues(2*num_eigen))
  allocate(eindex(2*num_eigen))

  X = TpetraMultiVector(map, num_eigen)

  ! Step 0: crate a handle
  eigen_handle = TrilinosEigenSolver(); FORTRILINOS_CHECK_IERR()

  ! ------------------------------------------------------------------
  ! Explicit setup and solve
  ! ------------------------------------------------------------------
  ! Step 1: initialize a handle
  call eigen_handle%init(comm); FORTRILINOS_CHECK_IERR()

  ! Step 2: setup the problem
  call eigen_handle%setup_matrix(A); FORTRILINOS_CHECK_IERR()
  call eigen_handle%setup_matrix_rhs(B); FORTRILINOS_CHECK_IERR()
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

#if FORTRILINOS_USE_MPI
  ! Finalize MPI must be called after releasing all handles
  call MPI_FINALIZE(ierr)
#endif


end program
